FOLDERS, IDS, RESTS, = glob_wildcards("{dir, (\w*\/)?}{id, [^_]+}_{rest, \w*_R2_\w*|\w*2}.fastq")
DATABASE = "Greengene-13-5-99.fasta"
CDHIT = "/home/manager/Folder_shared/cd-hit-v4.6.8-2017-0621"
INPUT_FOLDER = "/home/manager/test_reads"
TRIMMOMATIC = "trimmomatic-0.32.jar"
SIMILARITY_CUTOFF = 0.97
p = 170
q = 100
THREADS = 4
MINLEN = 80

rule all:
    input:
        "pooled/OTU.txt"


rule prepare:
    output:
        "sample1",
        "sample2",
        "SAMPLE_FILE",
        expand("{id}/{id}_1.fastq", id=IDS),
        expand("{id}/{id}_2.fastq", id=IDS)
    run:
        import os
        import random
        import re
        from shutil import copyfile

        def make_folder(filename, match_obj, input_dir, samples, postfix):
            sample_name = match_obj.group(1)
            sample_dir_name = '{}/{}/'.format(input_dir, sample_name)
            short_filename_with_postfix = sample_name + postfix
            if sample_name not in samples:
                samples[sample_name] = []
                if not os.path.exists(sample_dir_name):
                    os.makedirs(sample_dir_name)
            samples[sample_name].append(short_filename_with_postfix)
            os.rename('{}/{}'.format(input_dir, filename), sample_dir_name + short_filename_with_postfix)

        def put_sample_to_folder(filename, input_dir, samples):
            match1_obj = re.match(r"([^_]+)_(\w*_R1_\w*|\w*_R1|1).fastq", filename)
            if match1_obj:
                postfix = "_1.fastq"
                make_folder(filename, match1_obj, input_dir, samples, postfix)
            else:
                match2_obj = re.match(r"([^_]+)_(\w*_R2_\w*|\w*_R2|2).fastq", filename)
                if match2_obj:
                    postfix = "_2.fastq"
                    make_folder(filename, match2_obj, input_dir, samples, postfix)

        samples = {}
        input_dir = os.getcwd()
        for filename in os.listdir(input_dir):
            put_sample_to_folder(filename, input_dir, samples)

        sample, names = random.choice(list(samples.items()))
        names.sort()
        sample1 = '{}/{}/{}'.format(input_dir, sample, names[0])
        sample2 = '{}/{}/{}'.format(input_dir, sample, names[1])

        with open(input_dir + "/sample1", "w") as f:
            f.write(sample1)

        with open(input_dir + "/sample2", "w") as f:
            f.write(sample2)

        with open(input_dir + "/SAMPLE_FILE", "w") as f:
            for key, value in sorted(samples.items()):
                f.write(key + " " + " ".join(sorted(value)) + "\n")




rule spliced_db:
    input:
        sample1="sample1",
        sample2="sample2"
    output:
        "spliced_db-R1",
        "spliced_db-R2"
    shell:
        "{CDHIT}/usecases/Miseq-16S/16S-ref-db-PE-splice.pl -i `cat {input.sample1}` -j `cat {input.sample2}` -d {DATABASE} \
        -o spliced_db -p {p} -q {q} -c 0.99"


rule qc:
    input:
        "SAMPLE_FILE",
        R1="{id}/{id}_1.fastq",
        R2="{id}/{id}_2.fastq"        
    params:
        sample="{id}"
    output:
        R1_fq="{id}/qc/R1.fq",
        R2_fq="{id}/qc/R2.fq"
    shell:
        """
        java -jar {TRIMMOMATIC} PE -threads {THREADS} -phred33 {input.R1} {input.R2} {INPUT_FOLDER}/{params.sample}/qc/R1.fq \
        {INPUT_FOLDER}/{params.sample}/qc/R1-s.fq {INPUT_FOLDER}/{params.sample}/qc/R2.fq {INPUT_FOLDER}/{params.sample}/qc/R2-s.fq \
        SLIDINGWINDOW:4:20 LEADING:3 TRAILING:3 MINLEN:{MINLEN} MAXINFO:80:0.5 \
        1>{INPUT_FOLDER}/{params.sample}/qc/qc.stdout 2>{INPUT_FOLDER}/{params.sample}/qc/qc.stderr
        rm -f {params.sample}/qc/R1-s.fq {params.sample}/qc/R2-s.fq
        """


rule fasta_from_fastq:
    input:
        "{id}/qc/R1.fq",
        "{id}/qc/R2.fq"
    output:
        "{id}/qc/R1.fa",
        "{id}/qc/R2.fa"
    run:
        import os
        import random
        input_dir = os.getcwd()
        for ind, el in enumerate(output):
            count = 0
            with open(input_dir + "/" + input[ind]) as fq, open(input_dir + "/" + output[ind], "w") as fa:
                for line in fq:
                    count += 1
                    if count % 4 == 1:
                        fa.write(">Sample|" + input[0].split("/")[0] + "|" + str(count // 4 + 1) + " " + line)
                    elif count % 4 == 2:
                        fa.write(line)


rule seq_nr:
    input:
        R1_fa="{id}/qc/R1.fa",
        R2_fa="{id}/qc/R2.fa"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.nr",
        "{id}/otu/seq.nr.2",
        "{id}/otu/seq.nr.clstr"
    shell:
        """
        {CDHIT}/cd-hit-est -i {params.sample}/qc/R1.fa -j {params.sample}/qc/R2.fa -o {params.sample}/otu/seq.nr \
        -op {params.sample}/otu/seq.nr.2 -sf 1 -sc 1 -P 1 -r 0 -cx {p} -cy {q} -c 1.0  -n 10 -G 1 -b 1  -T 1 -M 8000 \
        -d 0 -p 1 > {params.sample}/otu/seq.nr.log
        """


rule chimeric_clstr:
    input:
        "{id}/otu/seq.nr",
        "{id}/otu/seq.nr.2"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.chimeric-clstr.R1",
        "{id}/otu/seq.chimeric-clstr.R2"
    shell:
        """
        {CDHIT}/cd-hit-est -i {params.sample}/otu/seq.nr   -o {params.sample}/otu/seq.chimeric-clstr.R1 -r 0 -cx 75 \
         -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > {params.sample}/otu/seq.chimeric-clstr.R1.log
        {CDHIT}/cd-hit-est -i {params.sample}/otu/seq.nr.2 -o {params.sample}/otu/seq.chimeric-clstr.R2 -r 0 -cx 75 \
        -c 0.99 -n 10 -G 0 -b 1 -A 50 -T 1 -M 8000  -d 0 -p 1 > {params.sample}/otu/seq.chimeric-clstr.R2.log
        """


rule seq_99:
    input:
        "{id}/otu/seq.nr",
        "{id}/otu/seq.nr.2"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.99",
        "{id}/otu/seq.99.2"
    shell:
        "{CDHIT}/cd-hit-est -i {params.sample}/otu/seq.nr -j {params.sample}/otu/seq.nr.2 -o {params.sample}/otu/seq.99 \
        -op {params.sample}/otu/seq.99.2 -P 1 -r 0 -cx {p} -cy {q} -c 0.99 -n 10 -G 1 -b 1  -T 1 -M 8000  -d 0 -p 1 \
        > {params.sample}/otu/seq.99.log"


rule seq_99f:
    input:
        "{id}/otu/seq.99",
        "{id}/otu/seq.99.2",
        "{id}/otu/seq.chimeric-clstr.R1",
        "{id}/otu/seq.chimeric-clstr.R2"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.99f",
        "{id}/otu/seq.99f.2"
    shell:
        "{CDHIT}/usecases/Miseq-16S/filter-chimeric-and-small.pl -c 0.0001 -k {params.sample}/otu/seq.nr.clstr \
        -i {params.sample}/otu/seq.chimeric-clstr.R1.clstr -j {params.sample}/otu/seq.chimeric-clstr.R2.clstr \
        -a {params.sample}/otu/seq.99.clstr -f {params.sample}/otu/seq.99 -g {params.sample}/otu/seq.99.2 -o {params.sample}/otu/seq.99f"


rule seq_99f_all:
    input:
        "{id}/otu/seq.nr",
        "{id}/otu/seq.99f"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.99f-all.clstr"
    shell:
        "{CDHIT}/clstr_rev.pl {params.sample}/otu/seq.nr.clstr {params.sample}/otu/seq.99f.clstr > {params.sample}/otu/seq.99f-all.clstr"

rule seq_97:
    input:
        "{id}/otu/seq.99f"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.97"
    shell:
        "{CDHIT}/cd-hit-est -i {params.sample}/otu/seq.99f -j {params.sample}/otu/seq.99f.2 -o {params.sample}/otu/seq.97 \
        -op {params.sample}/otu/seq.97.2 -P 1 -r 0 -cx {p} -cy {q} -c {SIMILARITY_CUTOFF} -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > \
        {params.sample}/otu/seq.97.log"


rule seq_97_ref:
    input:
        "{id}/otu/seq.97",
        db1="spliced_db-R1",
        db2="spliced_db-R2"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.97.ref"
    shell:
        "{CDHIT}/cd-hit-est-2d -i {params.sample}/otu/seq.97 -j {params.sample}/otu/seq.97.2 -i2 {INPUT_FOLDER}/{input.db1} \
        -j2 {INPUT_FOLDER}/{input.db2} -o {params.sample}/otu/seq.97.ref -op {params.sample}/otu/seq.97.ref.2 -P 1 -r 0 \
        -cx {p} -cy {q} -c {SIMILARITY_CUTOFF} -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > {params.sample}/otu/seq.97.ref.log"


rule seq_97_all:
    input:
        "{id}/otu/seq.99f-all.clstr",
        "{id}/otu/seq.97"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.97-all.clstr"
    shell:
        "{CDHIT}/clstr_rev.pl {params.sample}/otu/seq.99f-all.clstr {params.sample}/otu/seq.97.clstr > {params.sample}/otu/seq.97-all.clstr"


rule seq_97_reftop:
    input:
        "{id}/otu/seq.97.ref"
    params:
        sample="{id}"
    output:
        "{id}/otu/seq.97.reftop.clstr"
    shell:
        "{CDHIT}/usecases/Miseq-16S/filter-nontop-ref.pl < {params.sample}/otu/seq.97.ref.clstr > {params.sample}/otu/seq.97.reftop.clstr"


rule OTU_all:
    input:
        "{id}/otu/seq.97-all.clstr",
        "{id}/otu/seq.97.reftop.clstr"
    params:
        sample="{id}"
    output:
        "{id}/otu/OTU.clstr"
    shell:
        """
        {CDHIT}/clstr_merge.pl {params.sample}/otu/seq.97-all.clstr {params.sample}/otu/seq.97.reftop.clstr > {params.sample}/otu/OTU.clstr
        rm -f {params.sample}/otu/seq.chimeric-clstr.R1 {params.sample}/otu/seq.chimeric-clstr.R1.log \
        {params.sample}/otu/seq.chimeric-clstr.R2 {params.sample}/otu/seq.chimeric-clstr.R2.log
        rm -f {params.sample}/otu/seq.97.ref {params.sample}/otu/seq.97.ref.2 {params.sample}/otu/seq.97.ref.log
        mv {params.sample}/otu/seq.99f.log {params.sample}/otu/chimeric-small-clusters-list.txt
        """


rule pooling:
    input:
        expand("{id}/otu/OTU.clstr", id=IDS)
    output:
        "pooled/SAMPLE_FILE"
    shell:
        """
        {CDHIT}/usecases/Miseq-16S/pool_samples.pl -s SAMPLE_FILE -o pooled
        cp SAMPLE_FILE pooled
        """


rule pooled_seq_97:
    input:
        "pooled/SAMPLE_FILE"
    output:
        "pooled/seq.97"
    shell:
        """
        mkdir pooled/otu-pooled
        {CDHIT}/cd-hit-est -i pooled/seq.99f -j pooled/seq.99f.2 -o pooled/seq.97 -op pooled/seq.97.2 -P 1 -r 0 \
        -cx {p} -cy {q} -c {SIMILARITY_CUTOFF} -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 > pooled/seq.97.log
        """

rule pooled_seq_97_ref:
    input:
        "pooled/seq.97",
        db1="spliced_db-R1",
        db2="spliced_db-R2"
    output:
        "pooled/seq.97.ref"
    shell:
        "{CDHIT}/cd-hit-est-2d -i pooled/seq.97 -j pooled/seq.97.2 -i2 {INPUT_FOLDER}/{input.db1} -j2 {INPUT_FOLDER}/{input.db2} \
        -o pooled/seq.97.ref -op pooled/seq.97.ref.2 -P 1 -r 0 -cx {p} -cy {q} -c {SIMILARITY_CUTOFF} -n 10 -G 1 -b 10  -T 1 -M 8000  -d 0 -p 1 \
        > pooled/seq.97.ref.log"


rule pooled_seq_97_all:
    input:
        "pooled/seq.97",
    output:
        "pooled/seq.97-all.clstr"
    shell:
        "{CDHIT}/clstr_rev.pl pooled/seq.99f-all.clstr pooled/seq.97.clstr > pooled/seq.97-all.clstr"


rule pooled_seq_97_reftop:
    input:
        "pooled/seq.97.ref",
    output:
        "pooled/seq.97.reftop.clstr"
    shell:
        "{CDHIT}/usecases/Miseq-16S/filter-nontop-ref.pl < pooled/seq.97.ref.clstr > pooled/seq.97.reftop.clstr"


rule pooled_OTU_clstr:
    input:
        "pooled/seq.97-all.clstr",
        "pooled/seq.97.reftop.clstr"
    output:
        "pooled/OTU.clstr"
    shell:
        "{CDHIT}/clstr_merge.pl pooled/seq.97-all.clstr pooled/seq.97.reftop.clstr > pooled/OTU.clstr"


rule pooled_OTU_txt:
    input:
        "pooled/OTU.clstr"
    output:
        "pooled/OTU.txt"
    shell:
        "{CDHIT}/usecases/Miseq-16S/clstr_2_OTU_table.pl -i pooled/OTU.clstr -o pooled/OTU.txt"
