#! /usr/bin/perl

# Many thanks to Mayeul d'Avezac for providing this patch

move();

my @friends = ();
my @class = ();
my $template = "";
my $innertemplate = "";
my $name = "";
sub move()
{
  my $file = $_[0];

  my $isinclass=0;
  my $opened = 0;

  while( ($_ = <STDIN>)  )
  {
    if( /^(\s+|)class(\s+|)(\S+)/ and /\{(\s+)$/ )
    {
      /^(\s+|)class(\s+|)(\S+)/;
      $isinclass = 1;
      $name = $3;
      $name =~ s/://;
      push @class, ( $_ );
    }
    if( /template(\s+|)<((\S|\s)+)>/ and $isinclass == 0 )
    {
      $template = $2;
      $innertemplate = $template;
      while ( $innertemplate =~ s/class// ) {}
      while ( $innertemplate =~ s/types::t_int// ) {}
      while ( $innertemplate =~ s/int// ){}
      $innertemplate =~ s/^\s+//;
      $innertemplate =~ s/\s+$//;
    }
    $opened++ if( $isinclass and /{/ );
    if( /friend/ and not /friend(\s+)class/ )
      { grab_friend(); }
    elsif ( $isinclass == 1 ) 
      { push @class, ($_); }
    else { print $_; }
    if( $isinclass and /}/ )
    {
      $opened--;
      if(  $opened == 0 )
      {
        print_friends();
        $isinclass = 0;
        $template = "";
        $innertemplate = "";
        $name = "";
      }
    }

  }
}

sub grab_friend()
{
  return if ( /\)(\s+|);/ );
  my $line = $_;
  my $opened = 0;
  
  push @friends, ( sprintf "$1 template<%s>\n", $template  ) if( $template =~ /\S+/ );
  do
  {
    my $line = $_;

    while ( $line =~ s/$name(\s+|)<(.*)>/TEMPLATED$1<$2>/ ) {}
    if( $line =~ /$name/ )
    {
      $line =~ s/$name/$name<$innertemplate>/;
    }
    while ($line =~ s/TEMPLATED(\s+|)<(.+)>/$name$1<$2>/) {}

    push @friends, ( $line );

    ++$opened if( /{/ );
    --$opened if( /}/ );
    if( $opened == 0 )
    {
      push @friends, ( "\n" );
      return;
    }
  }
  while ( $_ = <STDIN> );
}

sub print_friends()
{
  if( $#friends == -1 )
  {
    for( my $i = 1; $i < scalar( @class ); $i++ )
    {
      print $class[$i];
    }
    $#friends = -1;
    $#class = -1;
    return;
  }
  # print forward class declaration
  print "  class $name; \n\n";
  # print forward declarations of friend function
  for( my $i = 0; $i < scalar( @friends ); $i++ )
  {
    $line = $friends[$i];
    next if( $line !~ /friend/ );

    my $dummy = $line;
    $dummy =~ /^(\s+)/;
    print "$1template<$template>\n"  if( $template =~ /\S/ );
    $line =~ s/friend//;
    $line =~ s/{//;
    $line =~ s/\)(\s+|)$/\);/;
    print $line, "\n";
  }

  # print start of class
  my $line = $class[0];
  my $dummy = $class[0];
  $dummy =~ /^(\s+)/;
  print "\n\n$1template<$template>\n" if( $template =~ /\S/ );
  print $class[0];
 
  # print friends
  for( my $i = 0; $i < scalar( @friends ); $i++ )
  {
    $line = $friends[$i];
    next if( $line !~ /friend/ );

    my $dummy = $line;
    $dummy =~ /^(\s+)/;
    $line =~ s/{//;
    $line =~ s/\(/ <> \(/;
    $line =~ s/\)(\s+|)$/\);/;
    print $line, "\n";
  }

  # print rest of class definition
  for( my $i = 2; $i < scalar( @class ); $i++ )
  {
    print $class[$i];
  }
  print "\n\n";
  # print friend function definition
  for( my $i = 0; $i < scalar( @friends ); $i++ )
  {
    $line = $friends[$i];
    $line =~ s/friend//;
    print $line;
  }
  $#friends = -1;
  $#class = -1;
}

