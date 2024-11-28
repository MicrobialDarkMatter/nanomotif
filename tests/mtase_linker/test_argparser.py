from nanomotif.argparser import create_parser

def test_argparse():
    parser = create_parser()

    
    args = parser.parse_args([
        "MTase-linker", "run",
        '--threads', '10', 
        '--forceall', 'True', 
        '--dryrun', 'True', 
        '--assembly', 'assembly_file', 
        '--contig_bin', 'contig_bin_file', 
        '--bin_motifs', 'motifs_file', 
        '-d', 'ML_dependencies_dir', 
        '--identity', '90', 
        '--qcovs', '90', 
        '--out', 'output_dir'
    ])

    # Check that the arguments are correctly parsed
    assert args.threads == 10
    assert args.forceall is True 
    assert args.dryrun is True
    assert args.assembly == 'assembly_file'
    assert args.contig_bin == 'contig_bin_file'
    assert args.bin_motifs == 'motifs_file'
    assert args.dependency_dir == 'ML_dependencies_dir'
    assert args.identity == '90' 
    assert args.qcovs == '90' 
    assert args.out == 'output_dir'
    
    print("All assertions passed!")

# Call the test function
test_argparse()