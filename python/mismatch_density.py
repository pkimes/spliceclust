import pysam

f = "/Users/pkimes/Downloads/sorted_region_chr6_test.bam"
samfile = pysam.Samfile(f, "rb")

print "###################"
# check different ways to iterate
print len(list(samfile.fetch()))
print len(list(samfile.fetch( "chr6", 10, 200 )))
print len(list(samfile.fetch()))

eg_pileups = samfile.pileup("chr6", 30000000, 40000000)
