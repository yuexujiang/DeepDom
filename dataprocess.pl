#!/usr/bin/perl
use Getopt::Long;
use Pod::Usage;
$seq_file = "";
$label_file = "";
$output_seq="";
$output_label = "";
$win=300;                        #specify the window length
$stride=80;                      #specify the stride length
$help=0;
$argnum = scalar(@ARGV);
GetOptions("input_seq=s" 		=> \$seq_file,
           "input_label=s" 		=> \$label_file,
           "output_seq=s"  		=> \$output_seq,
           "output_label=s"  		=> \$output_label,
           "windowsize=i" 		=> \$win,
           "stridesize=i" 		=> \$stride,
           "h|help" 		=> \$help) or pod2usage(-exitval => 2, -verbose => 2);

pod2usage(-exitval => 2,-verbose => 2) if ($help == 1);
pod2usage(-msg => "Invalid number of arguments!", -exitval => 2, -verbose => 2) if ($argnum < 4);
pod2usage(-msg => "please specify the name of your input sequence data(.fasta format)", -exitval => 2, -verbose => 2) if ($seq_file eq '');
pod2usage(-msg => "please specify the name of your input label data", -exitval => 2, -verbose => 2) if ($label_file eq '');
pod2usage(-msg => "please specify the output file name of the processed sequence data", -exitval => 2, -verbose => 2) if ($output_seq eq '');
pod2usage(-msg => "please specify the output file name of the processed label data", -exitval => 2, -verbose => 2) if ($output_label eq '');

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



__END__

=head1 NAME

dataprocess.pl

=head1 SYNOPSIS

perl dataprocess.pl [options] -input_seq <sequence file> -input_label <label file> -output_seq <output sequence file> -output_label <output label file>

=head1 ARGUMENTS

=over

=item B<-input_seq> <sequence file>

the name of your input sequence data, in the fasta format.

=item B<-input_label> <label file>

the name of your input label data.

=item B<-output_seq> <output sequence file>

the output file name of the processed sequence data

=item B<-output_label> <output label file>

the output file name of the processed label data

=back

=head1 OPTIONS

=over 
	   
=item B<-windowsize> <size of window> 

the size of overlapping window.(Default: 300)

=item B<-stridesize> <size of stride> 

the size of stride.(Default: 80)

=item B<-h|--help>   

Show help information.

=back

=head1 OUTPUT

This program will generate processed sequence data and label data.

=cut

