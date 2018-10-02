#!/usr/bin/perl

open $uni_seq_tr, "<","xxxxxxxxxxxxxxxxxxxxxxxxxx"or die "Can't open file: $!";                  #the positon of your sequence data (.fasta format)
                
open $cath_lable_tr, "<","xxxxxxxxxxxxxxxxxxxxxxx"or die "Can't open file: $!";                  #the positon of your label data

open $uni_seq_ntr, ">","xxxxxxxxxxxxxxxxxxxxxxxxx"or die "Can't open file: $!";         #the output position of the processed sequence data

open $cath_lable_ntr, ">","xxxxxxxxxxxxxxxxxxxxxx"or die "Can't open file: $!";         #the output position of the processed label data



$win=300;                        #specify the window length
$stride=80;                      #specify the stride length

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



