
import subprocess

def gtf2bed(file, gtf2bedscriptconvertor, utr=False):
    '''
    Function description:
        This function convert a gtf file into bed file with 12 columns (bed12)
        with utr if exists (.tmp) and without utr (.bed).
    
    Input:
        - file: a gtf file to convert in bed format.
        - utr: if you want to del utr information choose True.
    '''    
    path_to_script = "./"+gtf2bedscriptconvertor
    
    if utr == True:
        out = file+".bed"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
    
    if utr == False:
#        if "_reference" in file:
        out = file+".tmp"
#        else:
#            out = file+".bed"
        cmd = "{0} {1}  >> {2}".format(path_to_script, file, out)
        subprocess.call(cmd, shell=True)
#        if ".tmp" in out:
            
#            file_out = out+".bed"
#            path = "awk '{OFS='\t';split($11,a,','); split($12,b,','); A=''; B=''; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A''size',';B=B''start',';} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}'"
#            cmd = "{0} >> {1}".format(path, file_out)
#            subprocess.call(cmd, shell = True)
    #         cmd = "awk '{OFS='\t';split($11,a,','); split($12,b,','); A=''; B=''; if($7==$8) next; j=0; for(i=1;i<length(a);i++) if(($2+b[i]+a[i])>$7 && ($2+b[i])<$8) {j++; start=$2+b[i]-$7; size=a[i]; if(($2+b[i])<=$7) {start=0;size=size-($7-($2+b[i]));} if(($2+a[i]+b[i])>=$8) {size=size-($2+a[i]+b[i]-$8);} A=A''size',';B=B''start',';} print $1,$7,$8,$4,$5,$6,$7,$8,$9,j,A,B;}' out > file.bed"
        with open(out, 'r') as file_tmp:
            with open(file+".bed", 'w') as file_bed:
                for line in file_tmp.read().splitlines():
                    line = line.split('\t')
#                        print(line)
                    a = line[10].split(',')
#                        print(a)
                    b = line[11].split(',')
                    A = ""
                    B = ""
                    size = ""
                    start = ""
#                        
                    if line[6] == line[7]:
                        pass
                    j=0
                    for i in range(len(a)-1):
                        if (int(line[1])+int(b[i])+int(a[i])) > int(line[6]) and (int(line[1])+int(b[i])) < int(line[7]):
                            j+=1
                            start = int(line[1])+int(b[i])-int(line[6])
                            size = int(a[i])
                            if (int(line[1])+int(b[i]) <= int(line[6])):
                                start = 0
                                size = size-(int(line[6])-(int(line[1])+int(b[i])))
                        if (int(line[1])+int(a[i])+int(b[i])) >= int(line[7]):
                            size = size-(int(line[1])+int(a[i])+int(b[i])-int(line[7]))
                        A = A+''+str(size)+','
                        B = B+''+str(start)+','
                        if j == int(line[9]):
                            new_line = line[0]+'\t'+line[6]+'\t'+line[7]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6]+'\t'+line[7]+'\t'+line[8]+'\t'+str(j)+'\t'+str(A)+'\t'+str(B)
                            file_bed.write(new_line+'\n')
                            
                            
if __name__ == "__main__":
    gtf_file = "home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/predicted_transcripts_validation_with_new_annotations/remapped_mouse_cgalcode_tot.gtf"

    script = "/home/niguilla/Documents/Software/gtf_validator_of_cgalcode_prediction/dependencies/gtf2bed.pl"

    
    gtf2bed(gtf_file, script, False)