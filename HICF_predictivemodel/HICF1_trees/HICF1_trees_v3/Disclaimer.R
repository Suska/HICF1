#Disclaimer function

#This produces a pdf that can be glued to any preliminary data report

pdf("Disclaimer.pdf")
textplot("Analysis by Dr. Susanne Weller\n\nThese analyses is still preliminary!
\nGraphs should not be used for presentations and publications!\n", family="Helvetica")
dev.off()