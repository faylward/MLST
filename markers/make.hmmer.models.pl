$input = $ARGV[0];

opendir MYDIR, "$input";
@contents = readdir MYDIR;

$s = 0;
for $files (@contents) {

	if($files =~ /\.align\.fasta/) {
	$s++;
		system("hmmbuild --informat afa $input/$files\.hmm $input/$files");
	print "Processing $files\n";
	}

	else {

next;



	}
	}

closedir MYDIR;

