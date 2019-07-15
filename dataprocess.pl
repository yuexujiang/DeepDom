#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
$seq_file;
$label_file;
$output_seq;
$output_label;
$win=200;                        #specify the window length
$stride=80;                      #specify the stride length
$help=0;
$argnum = scalar(@ARGV);
GetOptions("input_seq=s" 		=> \$seq_file,
           "input_label=s" 		=> \$label_file,
           "output_seq=s"  		=> \$output_seq,
           "output_label=s"  		=> \$output_label,
           "windowsize=i" 		=> \$win,
           "stridesize=i" 		=> \$stride,
           "h|help" 		=> \$help) or die "$0 [--help] for help\n";
if($help)
{
my $message = <<'END_MESSAGE';
======================================

Note:
For predicting, please make sure the windowsize and the stridesize consistant with the trained model
For training, please specify the label files for both input and output 

Usage:
[--input_seq FILENAME]     specify the input sequence file (required)
[--input_label FILENAME]   specify the label file for input sequence (optional, only used for training a new model)
[--output_seq FILENAME]    specify the output file name of processed sequence  (required)
[--output_label FILENAME]  specify the output file name of processed label file (optional, only used for training a new model)
[--windowsize NUMBER]      specify the windowsize (required, default is 200)
[--stridesize NUMBER]      specify the stridesize (required, default is 80)
 
==========================================
END_MESSAGE
 
print $message;

}
else
{
if(defined $seq_file && defined $output_seq && defined $label_file && defined $output_label)
{
open $uni_seq_tr, "<".$seq_file or die "Can't open file: $!";                  #the positon of your sequence data (.fasta format)

open $cath_lable_tr, "<".$label_file or die "Can't open file: $!";                  #the positon of your label data

open $uni_seq_ntr, ">".$output_seq or die "Can't open file: $!";         #the output position of the processed sequence data

open $cath_lable_ntr, ">".$output_label or die "Can't open file: $!";         #the output position of the processed label data


$id="";
while($row=<$uni_seq_tr>)
{
	$row=~s/\r//g;
	chomp $row;	

	if($row=~/>(.*)$/)
	{
		$id=$1;
		next;
	}

	$l=length($row);
  if($l<$win)
  {
    $name=$id."_".$l."_0left".$l;
    print $uni_seq_ntr ">".$name."\n";
    print $uni_seq_ntr $row;
    $dis=$win-$l;
    for($i=0;$i<$dis;$i++)
    {
      print $uni_seq_ntr "-";
    }
    print $uni_seq_ntr "\n";
  }
  else
  {
    $iter=int(($l-$win)/$stride+1);
    $left=($l-$win)%$stride;
    for($i=0;$i<$iter;$i++)
    {
      $name=$id."_".$l."_".$i;
      print $uni_seq_ntr ">".$name."\n";
      print $uni_seq_ntr substr($row,$i*$stride,$win)."\n";
    }
    if($left!=0)
    {
      $comp=$l-$iter*$stride;
      
      $name=$id."_".$l."_".$iter."left".$comp;
      print $uni_seq_ntr ">".$name."\n";
      print $uni_seq_ntr substr($row,$iter*$stride,$comp);
      $dis=$win-$comp;
      for($i=0;$i<$dis;$i++)
      {
        print $uni_seq_ntr "-";
      }
      print $uni_seq_ntr "\n";
    }
  }
}

while($row=<$cath_lable_tr>)
{
	$row=~s/\r//g;
	chomp $row;	

	if($row=~/>(.*)$/)
	{
		$id=$1;
		next;
	}
 
	$l=length($row);
  if($l<$win)
  {
    $name=$id."_".$l."_0left".$l;
    print $cath_lable_ntr ">".$name."\n";
    print $cath_lable_ntr $row;
    $dis=$win-$l;
    for($i=0;$i<$dis;$i++)
    {
      print $cath_lable_ntr "-";
    }
    print $cath_lable_ntr "\n";
  }
  else
  {
    $iter=int(($l-$win)/$stride+1);
    $left=($l-$win)%$stride;
    for($i=0;$i<$iter;$i++)
    {
      $name=$id."_".$l."_".$i;
      print $cath_lable_ntr ">".$name."\n";
      print $cath_lable_ntr substr($row,$i*$stride,$win)."\n";
    }
    if($left!=0)
    {
      $comp=$l-$iter*$stride;
      $name=$id."_".$l."_".$iter."left".$comp;
      print $cath_lable_ntr ">".$name."\n";
      print $cath_lable_ntr substr($row,$iter*$stride,$comp);
      $dis=$win-$comp;
      for($i=0;$i<$dis;$i++)
      {
        print $cath_lable_ntr "-";
      }
      print $cath_lable_ntr "\n";
    }
  }

}

close $uni_seq_ntr;

close $cath_lable_ntr;

}
elsif(defined $seq_file && defined $output_seq)
{
open $uni_seq_tr, "<".$seq_file or die "Can't open file: $!";                  #the positon of your sequence data (.fasta format)
open $uni_seq_ntr, ">".$output_seq or die "Can't open file: $!";         #the output position of the processed sequence data
$id="";
while($row=<$uni_seq_tr>)
{
	$row=~s/\r//g;
	chomp $row;	

	if($row=~/>(.*)$/)
	{
		$id=$1;
		next;
	}

	$l=length($row);
  if($l<$win)
  {
    $name=$id."_".$l."_0left".$l;
    print $uni_seq_ntr ">".$name."\n";
    print $uni_seq_ntr $row;
    $dis=$win-$l;
    for($i=0;$i<$dis;$i++)
    {
      print $uni_seq_ntr "-";
    }
    print $uni_seq_ntr "\n";
  }
  else
  {
    $iter=int(($l-$win)/$stride+1);
    $left=($l-$win)%$stride;
    for($i=0;$i<$iter;$i++)
    {
      $name=$id."_".$l."_".$i;
      print $uni_seq_ntr ">".$name."\n";
      print $uni_seq_ntr substr($row,$i*$stride,$win)."\n";
    }
    if($left!=0)
    {
      $comp=$l-$iter*$stride;
      
      $name=$id."_".$l."_".$iter."left".$comp;
      print $uni_seq_ntr ">".$name."\n";
      print $uni_seq_ntr substr($row,$iter*$stride,$comp);
      $dis=$win-$comp;
      for($i=0;$i<$dis;$i++)
      {
        print $uni_seq_ntr "-";
      }
      print $uni_seq_ntr "\n";
    }
  }
}
close $uni_seq_ntr;
}
else
{
print "$0 [--help] for help\n";
}
}

