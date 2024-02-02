using LazyArtifacts


function process_QE_output(QE_output_string)
	QE_iterations = parse.(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/convergence has been achieved/s/.[^0-9]*\([0-9]\+\).*/\1/p'`), String),
			    '\n')[1:end-1]) # Remove last empty character.
	QE_n_kpoints = parse(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/number of k points=/s/.[^0-9]*\([0-9]\+\).*/\1/p'`), String),
			    '\n')[1])
	QE_optim_steps = parse(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/.*\([0-9]\+\) bfgs steps*/p'`), String),
			    '\n')[1])
	QE_symmetry_metadata = read(pipeline(
		`echo $QE_output_string`,
		`sed -n '/Sym\. Ops\./p'`), String)
	QE_n_G_vectors = parse.(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
	       `sed -n 's/.[^0-9]*\([0-9]\+\) G-vectors.*/\1/p'`), String),
			    '\n')[1])
	(; QE_iterations, QE_n_kpoints, QE_optim_steps, QE_symmetry_metadata, QE_n_G_vectors)
end
