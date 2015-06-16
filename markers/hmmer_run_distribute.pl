$input = $ARGV[0];
$output = $ARGV[1];
$db = $ARGV[2];


opendir MYDIR, "$input";
@contents = readdir MYDIR;

$s = 0;
for $files (@contents) {

	if($files =~ /\.faa/) {
	$s++;
		system("hmmsearch -E 1e-5 --cpu 3 --tblout $output/$files\_hmmersearch $db $input/$files");
	print "Processing $files\n";
	}

	else {

next;



	}
	}

closedir MYDIR;

