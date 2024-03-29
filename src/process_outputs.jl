using LazyArtifacts


function process_QE_output(QE_output_string; computation="scf")
	scf_iterations = parse.(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/convergence has been achieved/s/.[^0-9]*\([0-9]\+\).*/\1/p'`), String),
			    '\n')[1:end-1]) # Remove last empty character.
	n_kpoints = parse(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/number of k points=/s/.[^0-9]*\([0-9]\+\).*/\1/p'`), String),
			    '\n')[1])
	symmetry_metadata = read(pipeline(
		`echo $QE_output_string`,
		`sed -n '/Sym\. Ops\./p'`), String)
	n_symmetries = parse(Int64, read(pipeline(
		`echo $symmetry_metadata`,
		`sed -n 's/.[^0-9]*\([0-9]\+\) Sym\. Ops\..*/\1/p'`), String))
	n_G_vectors = parse.(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
	       `sed -n 's/.[^0-9]*\([0-9]\+\) G-vectors.*/\1/p'`), String),
			    '\n')[1])
	wall_time_str = read(pipeline(
		`echo $QE_output_string`,
		`sed -n 's/.*PWSCF.[^0-9]*\(.*\) CPU.[^0-9]*\(.*\) WALL$/\2/p'`), String)
	# Extract minutes and seconds.
	matched = match(r"(?:(\d+)h)?(?:(\d+)m)?(\d+\.\d+)s", wall_time_str)
	seconds = parse(Float64, matched.captures[3])
	if matched.captures[1] !== nothing
		hours = parse(Int, matched.captures[1])
	else hours = 0
	end
	if matched.captures[2] !== nothing
		minutes = parse(Int, matched.captures[2])
	else minutes = 0
	end
	wall_time = 60 * 60 * hours + 60 * minutes  + seconds

	if computation == "relax"
		optim_steps = parse(Int64, split(
			    read(
			        pipeline(
		`echo $QE_output_string`,
		`sed -n '/.*\([0-9]\+\) bfgs steps*/p'`), String),
			    '\n')[1])
	else optim_steps = nothing
	end
	(; scf_iterations, n_kpoints, optim_steps, n_symmetries, symmetry_metadata,
	 n_G_vectors, wall_time)
end

function process_DFTK_output(output)
	# We assume the timer has been reset properly before launch.
	scf_iterations = output.scfres.n_iter
	n_kpoints = length(output.scfres.basis.kpoints)
	n_G_vectors = length(output.scfres.basis.G_vectors)
	n_symmetries = length(output.scfres.basis.symmetries)
	
	(; scf_iterations, n_kpoints, n_symmetries, n_G_vectors, output.timings)
end
