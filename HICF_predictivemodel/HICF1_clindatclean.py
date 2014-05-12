#cleaning clinical data for HICF1 project

# Dr. Susanne Weller, 08/04/2014

filename = "clinicaldata_02052014.csv"
out = "clindat_R.txt"
output = open(out, "w")

print filename
with open(filename) as clindat:
        #next(clindat) #ommits header
        print clindat 
        for row in clindat:        
            print row
            #f = row.split("\t")
            #f[0].strip()
            #new_row = "\t".join(f)
            #output.write(new_row)
                       
    
output.close()  

