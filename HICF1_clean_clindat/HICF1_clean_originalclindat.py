# cleaning clinical data for HICF1 project
#This cleans the original clinical data file, especially the file names
#with unknown number of spaces behind and unpadded numbers

# Dr. Susanne Weller, 12/05/2014 (bypdw - day!)

filename = "HICF1_clinicaldata.txt"
out = "HICF1_clindat_clean11052014.txt"
output = open(out, "w")

with open(filename, 'rU') as clindat:
        #next(clindat) #ommits header 
        firstrow=True
        for row in clindat:        
            f = row.split("\t")
            if firstrow:
                f[0]="trialID"
                f[1]="trial"
                f[2]="dob"
                f[3]="gender"
                f[4]="dor"
                f[5]="aar"
                f[6]="vh"
                f[7]="vh3"
                f[8]="vh_mut"
                f[9]="13q14del"
                f[10]="freq13q14"
                f[11]="tri12"
                f[12]="freqtri"
                f[13]="11q23del"
                f[14]="freq11q23"
                f[15]="17pdel"
                f[16]="freq17pdel"
                f[17]="cd38"
                f[18]="mrd"
                f[19]="response"
                firstrow = False
            
            else:
                
                trialid = int(f[0])          
                    
                if  trialid < 10:
                    
                    f[0]= f[1]+"00" +f[0]
                    
                
                elif 10 <= trialid <= 100: 
                        f[0]=f[1]+"0" + f[0]
                       
                        
                else:
                    f[0]=f[1]+f[0]
                    
            
            new_row = "\t".join(f)
            output.write(new_row)
                       
    
output.close() 

