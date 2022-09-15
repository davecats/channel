use File::Copy;

open my $in,  '<', "dns.in"     or die "Can't read old file: $!";
open my $out, '>', "dns.in.new" or die "Can't write new file: $!";

while( <$in> )   # print the lines before the change
    {
        print $out $_;
        last if $. == 11; # line number before change
    }

my $line = <$in>;
$line = $ARGV[0]."\n";
print $out $line;

while( <$in> )   # print the rest of the lines
    {
        print $out $_;
    }

move("dns.in.new","dns.in");
