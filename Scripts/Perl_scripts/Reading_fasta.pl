open ($file,"C:/Users/HP/Downloads/Siddharth/Data/Chromosome_11_sample_data/chr11.fa") or die "Chromosome 11 file could not be opened";
$empty = '';
while (my $line = readline($file)){
  $empty += $line;
};
print "$empty is my collection of lines";