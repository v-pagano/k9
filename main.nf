nextflow.enable.dsl=2

@Grab('com.xlson.groovycsv:groovycsv:1.3')
import groovy.json.JsonBuilder
import groovy.json.JsonSlurper
import com.xlson.groovycsv.CsvParser

def readJsonFile(jsonFile) {
    def jsonSlurper = new JsonSlurper()
    def tupleFile = new File(jsonFile)
    String tupleString = tupleFile.text
    def tupleJson = jsonSlurper.parseText(tupleString)
    dictTSV = [:]
    sampleTSV = [:]

    def r1 = []
    def r2 = []
    def rg = []
    def nr = []
    def l = []

    for (i = 0; i < tupleJson['dataFiles'].size; i++) {
        if (tupleJson['dataFiles'][i]['fastqCode'] == 'R1') {
            r1.add(tupleJson['dataFiles'][i]['fastqPath'])
            nr.add(tupleJson['dataFiles'][i]['numberOfReads'])
            String strRet = '@RG\\t'
            strRet = strRet + 'ID:' + tupleJson['dataFiles'][i]["rgid"] + '\\t' 
            strRet = strRet + 'LB:' + tupleJson['dataFiles'][i]["rglb"] + '\\t'
            strRet = strRet + 'PU:' + tupleJson['dataFiles'][i]["rgpu"] + '\\t' 
            strRet = strRet + 'SM:' + tupleJson['dataFiles'][i]["rgsm"] + '\\t' 
            strRet = strRet + 'PL:' + tupleJson['dataFiles'][i]["rgpl"] + '\\t' 
            strRet = strRet + 'CN:tgen\\t'
            strRet = strRet + 'PM:' + tupleJson['dataFiles'][i]["rgpm"] + '\\t' 
            strRet = strRet + 'BC:Unknown'
            rg.add(strRet)
            m = [:]
            m = [subject: tupleJson['dataFiles'][i]['assayCode'], sex: tupleJson['sex'], sample: tupleJson['dataFiles'][i]['sampleName'], lane: tupleJson['dataFiles'][i]['rgpu'], fastq1: tupleJson['dataFiles'][i]['fastqPath'], rg: strRet]
            // Build out the tsv file for samples
            dictTSV[tupleJson['dataFiles'][i]['rgid']] = [:]
            dictTSV[tupleJson['dataFiles'][i]['rgid']] = [subject: tupleJson['dataFiles'][i]['assayCode'], sex: tupleJson['sex'], sample: tupleJson['dataFiles'][i]['sampleName'], lane: tupleJson['dataFiles'][i]['rgpu'], fastq1: tupleJson['dataFiles'][i]['fastqPath'], rg: strRet]
            sampleTSV[tupleJson['dataFiles'][i]['sampleName']] = []
            if (tupleJson['dataFiles'][i]['subGroup'] == 'Tumor') {
                dictTSV[tupleJson['dataFiles'][i]['rgid']]['status'] = 1
                m.status = 1
            } else {
                dictTSV[tupleJson['dataFiles'][i]['rgid']]['status'] = 0
                m.status = 0
            }

        }
        if (tupleJson['dataFiles'][i]['fastqCode'] == 'R2') {
            r2.add(tupleJson['dataFiles'][i]['fastqPath'])

            // Build out the tsv file for samples
            dictTSV[tupleJson['dataFiles'][i]['rgid']]['fastq2'] = tupleJson['dataFiles'][i]['fastqPath']
            m.fastq2 = tupleJson['dataFiles'][i]['fastqPath']
            l.add(m)
        }
    }

    outC = []
    for (i = 0; i < r1.size; i++) {
        outC.add(r1: r1[i], r2: r2[i], rg: rg[i], numReads: nr[i])
    }

    dictTSV.each { k, v -> sampleTSV[dictTSV[k]['sample']].add(v) }

    sampleList = []
    for (i in sampleTSV) {
        sampleList.add([sampleName: i.key, info: i.value])
    }
    outputTuple = [fileList: outC, sampleTSV: sampleList, rgTSV: dictTSV, fileList: l]
    return outputTuple
}

def typeOfInput(fileName) {
    def m = (fileName =~ /[^\.]*$/)
    if (m.size() > 0) {
        return m[0]
    } else {
        return ''
    }
}

def readCsvFile(fileName, fileType) {
    files = []
    sep = ','
    if (fileType == 'tsv') {
        sep = '\t'
    }
    data = new CsvParser().parse(
        new FileReader(fileName), 
        separator: ',', 
        quoteChar: "'", 
        readFirstLine: true,
        columnNames:['subject', 'sex', 'status', 'sample', 'lane', 'fastq1', 'fastq2', 'rg']
    )
    for (line in data) {
        files.add([
        subject: line.subject, 
        sex: line.sex, 
        status: line.status,
        sample: line.sample, 
        lane: line.lane, 
        fastq1: line.fastq1, 
        fastq2: line.fastq2,
        rg: line.rg
        ])
    }

    return files
}

include { MAPPING_TUMOR; MAPPING_NORMAL } from './workflows/mapping'
include { VARIANTCALLERS_NORMAL; VARIANTCALLERS_TUMOR } from './workflows/variantcallers'
include { ANNOTATION_NORMAL; ANNOTATION_TUMOR } from './workflows/annotation'
include { PB_HAPLOTYPECALLER; PB_DEEPVARIANT; PB_SOMATIC } from './workflows/parabricks'
include { PUBLISH_RESULTS } from './workflows/helpers'

workflow {
    inputType = typeOfInput(params.input)
    switch (inputType) {
        case 'json':
            fileList = readJsonFile(params.input).fileList
            break        
        case 'csv':
            fileList = readCsvFile(params.input, 'csv')
            break
        case 'tsv':
            fileList = readCsvFile(params.input, 'tsv')
            break
    }

    tumorList = fileList.findAll { it.status.toString() == '1' }
    normalList = fileList.findAll { it.status.toString() == '0' }
    tumorSample = tumorList.collect { it.sample }.unique()
    normalSample = normalList.collect { it.sample }.unique()

    publishList = Channel.empty()

    switch (params.pipeline_type) {
        case 'Alignment':
            if (tumorList.size() > 0) {
                MAPPING_TUMOR(tumorList)
                publishList = publishList.concat(MAPPING_TUMOR.out)
            }
            if (normalList.size() > 0) {
                MAPPING_NORMAL(normalList)
                publishList = publishList.concat(MAPPING_NORMAL.out)
            }
            break
        case 'Constitutional':
            if (tumorList.size() > 0) {
                MAPPING_TUMOR(tumorList)
                VARIANTCALLERS_TUMOR(MAPPING_TUMOR.out, tumorSample[0])
                ANNOTATION_TUMOR(VARIANTCALLERS_TUMOR.out, tumorSample[0])
                publishList = publishList.concat(MAPPING_TUMOR.out)
                publishList = publishList.concat(VARIANTCALLERS_TUMOR.out)
                publishList = publishList.concat(ANNOTATION_TUMOR.out)
            }
            if (normalList.size() > 0) {
                MAPPING_NORMAL(normalList)
                VARIANTCALLERS_NORMAL(MAPPING_NORMAL.out, normalSample[0])
                ANNOTATION_NORMAL(VARIANTCALLERS_NORMAL.out, normalSample[0])
                publishList = publishList.concat(MAPPING_NORMAL.out)
                publishList = publishList.concat(VARIANTCALLERS_NORMAL.out)
                publishList = publishList.concat(ANNOTATION_NORMAL.out)
            }
            break
        case 'Somatic':
            if (params.accelerated) {
                PB_SOMATIC(normalList, normalSample[0], tumorList, tumorSample[0])
                publishList = publishList.concat(PB_SOMATIC.out)
            }
            break

    }

    //Publish all of the files to the final output directory
    PUBLISH_RESULTS(publishList.flatten())
}