import os 
import gzip
import minus80 as m80

#------------------------------------------------------------------------------
#                "Labs grow great when old farts create workflows               
#                 whose ease they know they shall never enjoy."                 
#                                                           -Rob
#------------------------------------------------------------------------------

from src import utils
from pathlib import Path

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
FTP = FTPRemoteProvider()

from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
s3_key_id = os.environ.get('AWS_ACCESS_KEY')
s3_access_key = os.environ.get('AWS_SECRET_KEY')

S3 = S3RemoteProvider(
    endpoint_url='https://s3.msi.umn.edu',
    access_key_id=s3_key_id, 
    secret_access_key=s3_access_key
)

configfile: "config.yaml"

ncbi_id_map = {
    'NC_009144.3':'chr1',  'NC_009145.3':'chr2',  'NC_009146.3':'chr3',
    'NC_009147.3':'chr4',  'NC_009148.3':'chr5',  'NC_009149.3':'chr6',
    'NC_009150.3':'chr7',  'NC_009151.3':'chr8',  'NC_009152.3':'chr9',
    'NC_009153.3':'chr10', 'NC_009154.3':'chr11', 'NC_009155.3':'chr12',
    'NC_009156.3':'chr13', 'NC_009157.3':'chr14', 'NC_009158.3':'chr15',
    'NC_009159.3':'chr16', 'NC_009160.3':'chr17', 'NC_009161.3':'chr18',
    'NC_009162.3':'chr19', 'NC_009163.3':'chr20', 'NC_009164.3':'chr21',
    'NC_009165.3':'chr22', 'NC_009166.3':'chr23', 'NC_009167.3':'chr24',
    'NC_009168.3':'chr25', 'NC_009169.3':'chr26', 'NC_009170.3':'chr27',
    'NC_009171.3':'chr28', 'NC_009172.3':'chr29', 'NC_009173.3':'chr30',
    'NC_009174.3':'chr31', 'NC_009175.3':'chrX',  'NC_001640.1':'chrMt'
}

ensem_id_map = utils.gen_ensem_id_map()

# refseq - GCF v genbank - GCA
NCBI_ASSEM = config['NCBI']['ASSEMBLY']
ENSEMBL_ASSEM = config['ENSEMBL']['ASSEMBLY']
BUCKET = config['BUCKET']


# ----------------------------------------------------------
#       It's the all, ya'll
# ----------------------------------------------------------

rule all:
    input:
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_transcriptomic.nice.fna.gz')


# ----------------------------------------------------------
#       Make the LocusPocus Databases for the GFF/Fasta
# ----------------------------------------------------------

rule make_transcriptomic_fna:
    input:
        fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ncbiEquCab3/thawed/tinydb.json',
        gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ncbiEquCab3/thawed/tinydb.json'
    output:
        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_transcriptomic.nice.fna.gz',keep_local=True)
    run:
        import locuspocus as lp
        # Load the GFF and FNA dbs
        fna = lp.Fasta('ncbiEquCab3') 
        gff = lp.Loci('ncbiEquCab3')
        with open(output.fna,'w') as OUT:
            gff.set_primary_feature_type('gene')
            for gene in genes:
                longest = None
                max_length = 0
                # calulcate the length of each mRNA
                for feature in gene.subloci:
                    # skip non mRNA features
                    if feature.feature_type != 'mRNA':
                        continue
                    # calculate the total length of all the exons that make up the mRNA
                    exon_length = sum([len(x) for x in feature.subloci if x.feature_type == 'exon']) 
                    # Store info it its the longest
                    if exon_length > max_length:
                        longest = feature
                        max_length = exon_length
                if longest is None:
                    continue
                # Print out the nucleotides for the longest mRNA
                print(f">{gene.name}|{feature.name}",file=OUT)
                exon_seq = ''.join([fna[x.chromosome][x.start:x.end] for x in longest.subloci if x.feature_type == 'exon'])
                for chunk in [exon_seq[i:i+n] for i in range(0,len(exon_seq),90)]:
                    print(chunk,file=OUT)
                
                

rule make_locpoc_dbs:
    input:
        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz',keep_local=True),
        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz',keep_local=True)
    output:
        fna = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Loci.ncbiEquCab3/thawed/tinydb.json',
        gff = Path(m80.Config.cf.options.basedir) / 'datasets/v1/Fasta.ncbiEquCab3/thawed/tinydb.json'
    run:
        # Create the loci db
        import locuspocus as lp
        fna = lp.Fasta.from_file('ncbiEquCab3', input.fna)  
        # Create the GFF db
        gff = lp.Loci('ncbiEquCab3')
        gff.import_gff(input.gff)
        

# ----------------------------------------------------------
#       Make "nice" FASTAs
# ----------------------------------------------------------

rule nice_ncbi_fasta:
    input:
        fna = FTP.remote(f'{config["NCBI"]["FASTA"]}')
    output:
        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz',keep_local=True)
    run:
        with utils.RawFile(input.fna) as IN, open(output.fna,'w') as OUT:
            for line in IN:
                if line.startswith('>'):
                    name, *fields = line.lstrip('>').split()
                    if name in ncbi_id_map:
                        new_name = '>' + ncbi_id_map[name]
                        line = ' '.join([new_name, name] + fields + ['\n'])
                print(line,file=OUT,end='')


rule nice_ensembl_fasta:
    input:
        fna = FTP.remote(f'{config["ENSEMBL"]["FASTA"]}')
    output:
        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna.gz')
    run:
        with utils.RawFile(input.fna) as IN, open(output.fna,'w') as OUT:
            for line in IN:
                if line.startswith('>'):
                    name, *fields = line.lstrip('>').split()
                    if name in ensem_id_map:
                        new_name = '>' + ensem_id_map[name]
                        line = ' '.join([new_name, name] + fields + ['\n'])
                print(line,file=OUT,end='')


# ----------------------------------------------------------
#       Make "nice" GFFs
# ----------------------------------------------------------

rule nice_ncbi_gff:
    input:
        gff = FTP.remote(f'{config["NCBI"]["GFF"]}')
    output:
        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz',keep_local=True)
    run:
        with utils.RawFile(input.gff) as IN, \
            open(output.gff,'w') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in ncbi_id_map:
                    id = ncbi_id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')


rule nice_ensembl_gff:
    input:
        gff = FTP.remote(f'{config["ENSEMBL"]["GFF"]}')
    output:
        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff3.gz')
    run:
        with utils.RawFile(input.gff) as IN, \
            open(output.gff,'w') as OUT:
            for line in IN:
                id,*fields = line.split('\t')
                if id in ensem_id_map:
                    id = ensem_id_map[id]
                print(id,*fields,file=OUT,sep='\t',end='')


# ----------------------------------------------------------
#       Build STAR indices
# ----------------------------------------------------------

rule build_star_ncbi:
    input:
        gff = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.gff.gz'),
        fna = S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/{NCBI_ASSEM}_genomic.nice.fna.gz')
    output:
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/Genome'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/SA'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/SAindex'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrLength.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrName.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrNameLength.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/chrStart.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/exonGeTrInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/exonInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/geneInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/genomeParameters.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbInfo.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/sjdbList.out.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/transcriptInfo.tab')
    threads: int(f'{int(config["THREADS"]["STAR"])}')
    shell:
        f'''
          STAR \
          --runThreadN {{threads}} \
          --runMode genomeGenerate \
          --genomeDir {BUCKET}/public/refgen/{NCBI_ASSEM}/STAR_INDICES/ \
          --genomeFastaFiles {{input.fna}} \
          --sjdbGTFfile {{input.gff}} \
          --sjdbGTFtagExonParentTranscript Parent
        '''


rule build_star_ensembl:
    input:
        gff = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.gff3.gz'),
        fna = S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/{ENSEMBL_ASSEM}_genomic.nice.fna.gz')
    output:
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/Genome'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/SA'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/SAindex'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrLength.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrName.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrNameLength.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/chrStart.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/exonGeTrInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/exonInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/geneInfo.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/genomeParameters.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbInfo.txt'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbList.fromGTF.out.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/sjdbList.out.tab'),
        S3.remote(f'{BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/transcriptInfo.tab') 
    threads: int(f'{int(config["THREADS"]["STAR"])}')
    shell:
        f'''
          STAR \
          --runThreadN {{threads}} \
          --runMode genomeGenerate \
          --genomeDir {BUCKET}/public/refgen/{ENSEMBL_ASSEM}/STAR_INDICES/ \
          --genomeFastaFiles {{input.fna}} \
          --sjdbGTFfile {{input.gff}} \
          --sjdbGTFtagExonParentTranscript Parent 
        '''

# ----------------------------------------------------------
#       Build SALMON indices
# ----------------------------------------------------------

