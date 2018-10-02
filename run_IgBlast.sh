#1. Download IgBlast program and other required files
#IgBlast program can be downloaded from (https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST ).

#Other required files are also needed. Download the entire directory of internal_data and optional_file from (https://ftp.ncbi.nih.gov/blast/executables/igblast/release/ ).

#2. Make blast database for germline V, D, and J gene sequences.
#IgBlast allows you to search any germline databases of your choice (using -germline_db_V, -germline_db_J and -germline_db_D options).

#The NCBI mouse germline gene databases (i.e., mouse_gl_V, etc.) are supplied on our ftp site (see https://www.ncbi.nlm.nih.gov/igblast/ about database details).

#To search IMGT germline sequences, you need to download them from IMGT web site (http://www.imgt.org/vquest/refseqh.html#VQUEST ). You need to download all V, D and J sequences for whatever organisms you are interested in. Combine all V, all D and all J sequences, respectively, into separate files (i.e., one file for all V sequences, one for all D sequences and one file all for J sequences). After you have downloaded the sequences, invoke our utility tool edit_imgt_file.pl (download from the release/ directory) to process the sequences (to change the long IMGT definition lines to germline gene names only). 

#To search IMGT germline sequences, you need to download them from IMGT web site (http://www.imgt.org/vquest/refseqh.html#VQUEST ). You need to download all V, D and J sequences for whatever organisms you are interested in. Combine all V, all D and all J sequences, respectively, into separate files (i.e., one file for all V sequences, one for all D sequences and one file all for J sequences). After you have downloaded the sequences, invoke our utility tool edit_imgt_file.pl (download from the release/ directory) to process the sequences (to change the long IMGT definition lines to germline gene names only). For example:
#./edit_imgt_file.pl imgt_file > my_seq_file
#Then you can use NCBI’s makeblastdb tool to make the blast database from the output file. For example:

##Important links
##https://ncbi.github.io/igblast/cook/How-to-set-up.html
##https://changeo.readthedocs.io/en/version-0.3.7---igblast-1.7-fix/standard.html


###
# V-segment database
#perl edit_imgt_file.pl /Users/Pinedasans/VDJ_V2/IgBlast/Homo_sapiens/IG/IGHV.fasta > human_igh_v
#makeblastdb -parse_seqids -dbtype nucl -in human_igh_v
# D-segment database
#perl edit_imgt_file.pl /Users/Pinedasans/VDJ_V2/IgBlast/Homo_sapiens/IG/IGHD.fasta > human_igh_d
#makeblastdb -parse_seqids -dbtype nucl -in human_igh_d
# J-segment database
#perl edit_imgt_file.pl /Users/Pinedasans/VDJ_V2/IgBlast/Homo_sapiens/IG/IGHJ.fasta > human_igh_j
#makeblastdb -parse_seqids -dbtype nucl -in human_igh_j

FILES=/Users/Pinedasans/VDJ_V2/Data/FASTA/*.fasta
for f in $FILES
do
    echo $f
    var=$(basename $f)
    igblastn -germline_db_V human_igh_v -num_alignments_V 1 \
         -germline_db_D human_igh_d \
         -germline_db_J human_igh_j  -num_alignments_J 1 -auxiliary_data optional_file/human_gl.aux \
         -domain_system imgt -ig_seqtype Ig -organism human \
         -outfmt '7 std qseq sseq btop' \
         -query $f \
         -out $var.fmt7
         -num_alignments_V 1
        
    MakeDb.py igblast -i $var.fmt7 -s $f -r Homo_sapiens/IG/IGHV.fasta Homo_sapiens/IG/IGHD.fasta Homo_sapiens/IG/IGHJ.fasta --scores --cdr3 --regions
done