nextflow run https://github.com/v-pagano/k9 \
		 -name ecstatic_hamilton \
		 -params-file nf-1Xq2lSJgvuy3VQ.params.json \
		 -with-tower http://pnap-tower.tgen.org:8000/api \
		 -r master \
		 -profile dback \
		 -latest