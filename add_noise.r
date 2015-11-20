# add_noise.r                                       8/07 jlong@jimlong.org

# given some synthetic microarray data from a COPASI run, add simple noise.
# set the following variable to change the maximum deviation allowed, i.e.
# CLAMP = 0.25 means that the noise will not change a data point more than 
# 1/4 of its value 

CLAMP = 0.25

# this code invoked something like:
# R --vanilla --slave --args infile=<input_file> < add_noise.r > output_file
# where <input_file> has a commented-out header line followed by data, 
# and the first column is time, i.e:

# # Time     P7       P176     P198    P159    P191      P0      P1	
# 0          1        1	       1       1       1         1       1	
# 1          0.982247 2.63574  5.5281  2.28517 0.998523  3.5033  1.14449
# 2          0.966031 4.19939  9.66453 3.491   0.997198  5.88535 1.28263
# 3          0.951769 5.69412 13.4431  4.62131 0.99601   8.15203 1.41468
# 4          0.940614 7.12298 16.8948  5.6791  0.994945	10.3089  1.54091
# 5          0.933674 8.48887 20.048   6.66654 0.99399  12.3614  1.66159
# 6          0.931698 9.79455 22.9283  7.58491 0.993134 14.3144  1.77695	
# etc

# Example command: 
# R --vanilla --slave --args infile=500out.txt < add_noise.r > noisy_data.out
# Caution: takes a long time to run on large datasets

# read in data
flag=FALSE
args=unlist(strsplit(commandArgs(), "infile=", fixed=TRUE))
for(i in 1:length(args))
{
  if(args[i] == "")
  {
    file = args[i+1]
    flag = TRUE
    break
  }
}

if(flag)
{
  x = read.table(file, header=F)
}else
{
  cat("input file not found, exiting...\n")
  q()
}

num_of_genes = length(x[1,])  # number of cols in data
num_of_steps = length(x[,1])  # number of rows in data

# go through the matrix and add random noise
for(i in 1:num_of_steps)
{
  for(j in 2:num_of_genes)
  {
    # the amount of noise to add is picked from a normal distribution where
    # the mean is zero, and the standard deviation is 0.1*x[i,j]
    
    noise = rnorm(1, sd=0.1*x[i,j])
    
    # clamp noise at CLAMP of the value of x[i,j]
    while(abs(noise) > CLAMP * x[i,j])
      noise = 0.95 * noise
      
    cat(x[i,j] + noise,"")
  }
  cat("\n")
}

