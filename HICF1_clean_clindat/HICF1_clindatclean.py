#cleaning clinical data for HICF1 project

# Dr. Susanne Weller, 08/04/2014

filename = "clinicaldata_02052014.csv"
out = "clindat_R.txt"
output = open(out, "w")

with open(filename) as clindat:
        #next(clindat) #ommits header
        for row in clindat:        
                f = row.split("\t")
    
                print f[:3]
                f[0] = f[0].strip()
                print f[:3]
                
                new_row = "\t".join(f)
                output.write(new_row)
                       
    
output.close()  
