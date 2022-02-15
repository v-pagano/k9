include { publishResults } from '../process/helpers'

workflow PUBLISH_RESULTS {
    take: 
        f

    main:
        publishResults(f)

    emit:
        publishResults.out
}