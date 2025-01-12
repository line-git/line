#include <getopt.h>
#include <iostream>
#include <cstring>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <cmath>
// #include <complex>
#include "mpc.h"

using namespace std;

#include "global_vars.h"
#include "utils.h"
#include "setup.h"
#include "tensor_utils.h"
// #include "conversions.h"
#include "system_analyzer.h"
#include "eq_solver.h"
// #include "normalize.h"
#include "malloc_defs.h"
#include "poly_frac.h"
// #include "pf_normalize.h"
#include "codec.h"
#include "bc_gen.h"
#include "topology.h"
#include "kira.h"
extern "C" {
  // #include "jordan.h"
  // #include "algebra.h"
  #include "cpoly.h"
  #include "in_out.h"
  // #include "mp_roots_poly.h"
}

// global variables
int workprec;
int wp2;
int dbg = 0;
mpfr_t mpfr_tol;
double sleep_time = 0;
// double sleep_time = 0.2e6;
// double sleep_time = 0.01e6;
// double sleep_time = 0.005e6;


int main(int argc, char *argv[])
{
  int print = 0;

  //////
  // OPTIONS FROM TERMINAL
  //////
	int opt;
	int long_index = 0;

	// variables for options
	char *filepath_card = NULL;
	char *dir_parent = strdup("app");
	char *dir_bound = NULL;
	char *dir_run = NULL;
	char *filepath_result = NULL;
	int opt_bar = 0;
	int opt_one_eps = 0;
	int opt_eps_list = 0;
	int opt_checkpoint = 0;
	int opt_write = -3;
	int opt_kira_parallel = 1;
  int opt_kira_redo = -1;
  double opt_incr_prec = 1;
  int opt_prune_eps_abs = -20, opt_prune_eps_abs_from_terminal = 0;
  int opt_prune_eps_mode = 1, opt_prune_eps_mode_from_terminal = 0;

	// defining long options
	struct option long_options[] = {
		{"input", required_argument, NULL, 'i'},
		{"parent-dir", required_argument, NULL, 'p'},
		{"bound-dir", required_argument, NULL, 0},
		{"run-dir", required_argument, NULL, 0},
  	{"result-file", required_argument, NULL, 'r'},
  	{"bar", required_argument, NULL, 'b'},
		{"one-eps", no_argument, NULL, 0},
		{"eps-list", no_argument, NULL, 0},
		{"write", required_argument, NULL, 'w'},
		{"checkpoint", required_argument, NULL, 'c'},
		{"kira-parallel", required_argument, NULL, 0},
		{"kira-redo", required_argument, NULL, 0},
		{"incr-prec", required_argument, NULL, 0},
		{"prune-eps-abs", required_argument, NULL, 0},
		{"prune-eps-mode", required_argument, NULL, 0},
		{0, 0, 0, 0} // terminator
	};

	while ((opt = getopt_long(argc, argv, "i:p:r:bw:c:", long_options, &long_index)) != -1) {
		switch (opt) {
			case 'i':
				printf("Option --input (-i) has arg: %s\n", optarg);
				filepath_card = strdup(optarg);
				break;
			case 'p':
				printf("Option --parent-dir (-p) has arg: %s\n", optarg);
        free(dir_parent);
				dir_parent = strdup(optarg);
        break;
			case 'r':
				printf("Option --result-file (-r) has arg: %s\n", optarg);
				filepath_result = strdup(optarg);
				break;
			case 'b':
				printf("Option --bar (-b) activated.\n");
				opt_bar = 1;
				break;
			case 'w':
				printf("Option --write (-w) has arg: %s\n", optarg);
				opt_write = atoi(optarg);
        if (opt_write < -3 || opt_write > 2) {
          fprintf(stderr, "Invalid value for --write: %s. Must be integer in [-3, 2].\n", optarg);
          return 1;
        }
				break;
			case 'c':
				printf("Option --checkpoint (-c) has arg: %s\n", optarg);
				opt_checkpoint = atoi(optarg);
        if (opt_checkpoint < -1 || opt_checkpoint > 2) {
          fprintf(stderr, "Invalid value for --checkpoint: %s. Must be integer in [-1, 2].\n", optarg);
          return 1;
        }
				break;
			case 0: // long options without a short version
				if (strcmp("bound-dir", long_options[long_index].name) == 0) {
					printf("Option --bound-dir has arg: %s\n", optarg);
					dir_bound = strdup(optarg);
				} else if (strcmp("run-dir", long_options[long_index].name) == 0) {
					printf("Option --run-dir has arg: %s\n", optarg);
					dir_run = strdup(optarg);
				} else if (strcmp("one-eps", long_options[long_index].name) == 0) {
					printf("Option --one-eps activated.\n");
					opt_one_eps = 1;
				} else if (strcmp("eps-list", long_options[long_index].name) == 0) {
					printf("Option --eps-list activated.\n");
					opt_eps_list = 1;
				} else if (strcmp("kira-parallel", long_options[long_index].name) == 0) {
					char *endptr;
					long value = strtol(optarg, &endptr, 10);
					if (*endptr != '\0' || value <= 0) {
						fprintf(stderr, "Invalid value for --kira-parallel: %s. Must be a positive integer.\n", optarg);
						return 1;
					}
					printf("Option --kira-parallel has arg: %ld\n", value);
					opt_kira_parallel = (int)value;
				} else if (strcmp("kira-redo", long_options[long_index].name) == 0) {
					char *endptr;
					long value = strtol(optarg, &endptr, 10);
					if (*endptr != '\0' || !(value == 0 ||  value == 1 ||  value == -1)) {
						fprintf(stderr, "Invalid value for --kira-redo: %s. Must be integer in [-1, 1].\n", optarg);
						return 1;
					}
					printf("Option --kira-redo has arg: %ld\n", value);
					opt_kira_redo = (int)value;
        } else if (strcmp("incr-prec", long_options[long_index].name) == 0) {
					char *endptr;
					opt_incr_prec = strtod(optarg, &endptr);
					if (*endptr != '\0' || opt_incr_prec < 0 ) {
						fprintf(stderr, "Invalid value for --incr-prec: %s. Must be a positive double.\n", optarg);
						return 1;
					}
					printf("Option --incr-prec has arg: %f\n", opt_incr_prec);
        } else if (strcmp("prune-eps-abs", long_options[long_index].name) == 0) {
					char *endptr;
					long value = strtol(optarg, &endptr, 10);
					if (*endptr != '\0') {
						fprintf(stderr, "Invalid value for --prune-eps-abs: %s. Must be integer.\n", optarg);
						return 1;
					}
					printf("Option --prune-eps-abs has arg: %ld\n", value);
					opt_prune_eps_abs = (int)value;
          opt_prune_eps_abs_from_terminal = 1;
        } else if (strcmp("prune-eps-mode", long_options[long_index].name) == 0) {
					char *endptr;
					long value = strtol(optarg, &endptr, 10);
					if (*endptr != '\0' || !(value == 0 || value == 1 || value == -1)) {
						fprintf(stderr, "Invalid value for --prune-eps-mode: %s. Must be integer in [-1, 1].\n", optarg);
						return 1;
					}
					printf("Option --prune-eps-mode has arg: %ld\n", value);
					opt_prune_eps_mode = (int)value;
          opt_prune_eps_mode_from_terminal = 1;
        }
				break;
			case '?':
				printf("Unknown option: %s\n", argv[optind - 1]);
				return 1;
			case ':':
				printf("Missing argument for %s\n", argv[optind - 1]);
				return 1;
		}
	}

  //////
  // PROCESSING OPTIONS
  //////
  if (opt_checkpoint >= 1) {
    if (opt_kira_redo != 0) {
      cout << endl;
      cout << "Option --kira-redo = " << opt_kira_redo;
      cout << " while --checkpoint = " << opt_checkpoint << ";" << endl;
      cout << "setting --kira-redo to 0" << endl;
      opt_kira_redo = 0;
    }
  }

  // open virtual terminal for progress bars
  FILE *terminal = NULL;
  if (opt_bar) {
    // terminal = fopen("/dev/tty", "w");
    terminal = stderr;
  } else {
    terminal = fopen("/dev/null", "w");
    if (terminal == NULL) {
      perror("error while opening /dev/null.");
      exit(1);
    }
  }

  int opt_all_eps = 1;
  if (opt_one_eps) {
    opt_all_eps = 0;
  }

  //////
  // PARENT FILE PATHS
  //////
  cout << endl; cout << "open file card: " << filepath_card << endl;
  FILE *card_fptr = fopen(filepath_card, "r");
  if (card_fptr == NULL) {
    perror("error while opening file");
    exit(1);
  }
  char *dir_work = NULL;
  if (!read_param_str(&dir_work, (char*)"work-dir", card_fptr)) {
    perror("working directory not specified");
    fclose(card_fptr);
    exit(1);
  }
  fclose(card_fptr);
  if (dir_parent[strlen(dir_parent) -1] != '/') join_path(&dir_parent, dir_parent, (char*)"/");
  join_path(&dir_parent, dir_parent, dir_work);
  if (dir_parent[strlen(dir_parent) -1] != '/') join_path(&dir_parent, dir_parent, (char*)"/");

  char file_ext[] = ".txt";
  char *tmp_filepath = (char*) malloc(MAX_PATH_LEN*sizeof(char));
  char *tmp_filepath1 = (char*) malloc(MAX_PATH_LEN*sizeof(char));
  char *tmp_filepath2 = (char*) malloc(MAX_PATH_LEN*sizeof(char));
	char *dir_common = NULL;
	join_path(&dir_common, dir_parent, (char*)"common/");
	cout << "parent folder: " << dir_parent << endl;
	cout << "common dir: " << dir_common << endl;

  char *filepath_vars = NULL, *filepath_mats = NULL;
  char *filepath_branchcuts = NULL, *filepath_start = NULL;
  char *filepath_bound_behav = NULL, *filepath_bound_build = NULL;
  join_path(&filepath_vars, dir_common, (char*)"vars.txt");
  join_path(&filepath_branchcuts, dir_common, (char*)"branch_cuts.txt");
  join_path(&filepath_start, dir_common, (char*)"initial_point.txt");

  //////
  // LOAD SYMBOLS
  //////
  int ninvs;
	char **symbols;
	load_symbols(&ninvs, &symbols, filepath_vars);
	cout << "symbols: ";
	for (int s=0; s<=ninvs; s++) {
		printf("%s, ", symbols[s]);
	}
	printf("\n");

  // DETECT MASSES
  int linearize_mass;
  card_fptr = fopen(filepath_card, "r");
  if (!read_param_int(&linearize_mass, (char*)"lin-mass", card_fptr)) {
    cout << "lin-mass parameter not specified. Setting to defualt: 1." << endl;
    linearize_mass = 1;
  }

  int nmass, *mass, *is_mass;
  detect_masses(
    &nmass, &mass, &is_mass,
    ninvs, symbols
  );

  if (!linearize_mass) {
    nmass = 0;
    for (int s=0; s<ninvs; s++) {
      is_mass[s] = 0;
    }
  }

  cout << "num. kinematic invariants: " << ninvs << endl;
  cout << "num. masses: " << nmass << endl;
  cout << "invariants: ";
  for (int i=1; i<=ninvs; i++) {
    cout << symbols[i] << ", ";
  }
	cout << endl;
  cout << "masses: ";
  for (int m=0; m<nmass; m++) {
    cout << mass[m] << ", ";
  }
	cout << endl;
  cout << "is a mass?: ";
  for (int s=0; s<ninvs; s++) {
    cout << is_mass[s] << ", ";
  }
	cout << endl;

  //////
  // LOAD BRANCH CUTS
  //////
  // read branch cuts from file
  int nbranches;
  char **branch_cut_str;
  cout << endl; cout << "reading " << filepath_branchcuts << endl;
  load_branch_expr_str(&nbranches, &branch_cut_str, linearize_mass, filepath_branchcuts);
  cout << "BRANCH CUTS:" << endl;
  for (int b=0; b<nbranches; b++) {
    cout << branch_cut_str[b] << endl;
  }

  //////
  // LOAD PARAMETERS
  //////

  // EPSILON ABSOLUTE PRUNE
  cout << endl;
  if (opt_prune_eps_abs_from_terminal == 0) {
    fseek(card_fptr, 0, SEEK_SET);
    if (!read_param_int(&opt_prune_eps_abs, (char*)"prune-eps-abs", card_fptr)) {
      cout << "--prune-eps-abs option not specified. Using default: " << opt_prune_eps_abs << endl;
    } else {
      cout << "loaded --prune-eps-abs = " << opt_prune_eps_abs << " from card.";
      opt_prune_eps_abs_from_terminal = -1;
    }
  }
  int prune_eps_exp2 = opt_prune_eps_abs * log2(10);

  // EPSILON RELATIVE PRUNE
  if (opt_prune_eps_mode_from_terminal == 0) {
    fseek(card_fptr, 0, SEEK_SET);
    if (!read_param_int(&opt_prune_eps_mode, (char*)"prune-eps-mode", card_fptr)) {
      cout << "--prune-eps-mode option not specified. Using default: " << opt_prune_eps_mode << endl;
    } else {
      cout << "loaded --prune-eps-mode = " << opt_prune_eps_mode << " from card.";
      opt_prune_eps_mode_from_terminal = -1;
    }
  }

  if (opt_prune_eps_mode != -1) {
    if (
      opt_prune_eps_abs_from_terminal == 1 || \
      (opt_prune_eps_abs_from_terminal == -1 && opt_prune_eps_mode_from_terminal != 1)
    ) {
      perror("--prune-eps-abs specified while --prune-eps-mode != -1. Please omit --prune-eps-abs or choose --prune-eps-mode = -1");
      exit(1);
    }
  }

  // NLOOPS, ORDER, PRECISION, KINEMATICS

  // int eta_ord, eps_num, order = 5, precision = 16, nloops = 2;
  int eta_ord, eps_num;
  int nloops, order, precision;
  char *prev_run;
	// char **symbols;
	// int ninvs;
	char ***kin = new char**[2];
  load_parameters(
		&nloops, &order, &precision, &prev_run,
		&kin[1], &ninvs, &symbols,
		filepath_card
	);
  
  cout << endl; cout << "INPUT PARAMETERS:" << endl;
  cout << "nloops = " << nloops << endl;
  cout << "order = " << order << endl;
  cout << "precision = " << precision << endl;
  // if (prev_run) {
  //   cout << endl; cout << "prev_run: " << prev_run << endl;
  //   prepare_prev_run(&filepath_prev_run_PS, &filepath_bound, dir_parent, prev_run);
  //   cout << filepath_prev_run_PS << endl;
  //   cout << filepath_bound << endl;
  // }

	printf("target phase-space point: [\n");
	for (int s=0; s<ninvs; s++) {
		cout << "  " << symbols[s+1] << " = " << kin[1][s] << "," << endl;
	}
	printf("]\n");

	//////
	// TARGET HASH
	//////
	cout << endl; cout << "build target string..." << endl;
  char *target_str;
  char **param_keys = new char*[2];
  param_keys[0] = (char*)"ord";
  param_keys[1] = (char*)"prec";
  // [][20] = {"order", "precision"};
  int param_values[] = {order, precision};
	target_string(
		&target_str,
		2, (char**)param_keys, param_values,
		ninvs, symbols+1, kin[1]
	);
	cout << "target string: " << target_str << endl;

	cout << "generate hash..." << endl;
	char *target_hash;
	calculate_hash(&target_hash, target_str);
	cout << "target hash: " << target_hash << endl;

	//////
	// TARGET DIRECTORY
	//////
  char *dir_target = NULL;
	join_path(&dir_target, dir_parent, (char*)"points/");
  // create points directory
	if (make_dir(dir_target)) {
		printf("points directory already exists.\n");
	}

	join_path(&dir_target, dir_target, target_hash);
	join_path(&dir_target, dir_target, (char*)"/");
	cout << endl; cout << "target directory: " << dir_target << endl;

	// create target directory
	if (make_dir(dir_target)) {
		printf("hash target directory already exists.\n");
	}

  //////
  // SETUP
  //////
  char **eps_str;
  compute_setup(
    &workprec, &eta_ord, &eps_num, &eps_str,
    order, precision, nloops
  );
  wp2 = (10*workprec)/3; // binary working precision
  wp2 *= opt_incr_prec;

  mpfr_tol_set();

  cout << endl; cout << "INTERNAL PARAMETERS:" << endl;
  cout << "eta_ord = " << eta_ord << endl;
  cout << "workprec = " << workprec << endl;
  cout << "binary workprec = " << wp2 << endl;
  cout << "mpfr tolerance = "; mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
  
  cout << endl; cout << "EPSILON VALUES: " << endl;
  for (int ep=0 ; ep<eps_num ; ep++) {
    cout << eps_str[ep] << endl;
  }

  if (opt_eps_list) {
    exit(0);
  }

  if (opt_all_eps == 0) {
    eps_num = 1;
  }

  //////
  // CHOOSE STARTING POINT
  //////
  kin[0] = new char*[ninvs];
  // card_fptr = fopen(filepath_card, "r");

  // look for exit-sing option
  cout << endl;
  int exit_sing;
  fseek(card_fptr, 0, SEEK_SET);
  if (!read_param_int(&exit_sing, (char*)"exit-sing", card_fptr)) {
    cout << "exit-sing option not specified" << endl;
    exit_sing = 0;
  }
  cout << "exit-sing: " << exit_sing << endl;

  char *dir_amflow = NULL;
  int *MI_idx, dim_eta_less, dim;
  LI *MI_eta = NULL;
  if (exit_sing == -1) {
    //////
    // KIRA INTERFACE
    //////
    join_path(&dir_amflow, dir_parent, (char*)"amflow/");
    make_dir(dir_amflow);
    join_path(&dir_amflow, dir_amflow, target_hash);
    join_path(&dir_amflow, dir_amflow, (char*)"/");
    make_dir(dir_amflow);
    join_path(&dir_common, dir_amflow, (char*)"common/");

    call_kira(
      &MI_eta, &dim, &MI_idx, &dim_eta_less,
      opt_kira_redo, opt_kira_parallel, kin[1], dir_parent, dir_amflow,
      terminal
    );
  }
  // path to common files
  join_path(&filepath_mats, dir_common, (char*)"");
  join_path(&filepath_bound_behav, dir_common, (char*)"bound_behav.txt");
  join_path(&filepath_bound_build, dir_common, (char*)"bound_build.txt");
  
  // look for starting-point option
  int prev_run_start_pt = 0;
  fseek(card_fptr, 0, SEEK_SET);
  if (!read_PS_point(
    &kin[0],
    &ninvs, symbols+1, (char*)"starting-point", card_fptr
  )) {
    if (exit_sing == 0) {
      // #2BD: automatic selection of starting point
      perror(
        "neither exit-sing nor starting-point specified. Automatic selection of starting point not yet implemented"
      );
      fclose(card_fptr);
      exit(1);
    } else if (exit_sing == 1) {
      // choose to exit from singularity
      FILE *fptr = fopen(filepath_start, "r");
      if (!read_PS_point(
        &kin[0],
        &ninvs, symbols+1, (char*)"point", fptr
      )) {
        perror("error while reading file: PS point not found");
        fclose(fptr);
        exit(1);
      }
      fclose(fptr);
    }
  
  } else {
    if (exit_sing == 1) {
      cout << "both exit-sing = 0 and starting-point options specified. Overwrting default starting point." << endl;
      prev_run_start_pt = 0;
    } else if (exit_sing == 0) {
      prev_run_start_pt = 1;
    } else if (exit_sing == -1) {
      perror("both exit-sing = -1 and starting point options specified");
      exit(1);
    }

    printf("starting phase-space point: [\n");
    for (int s=0; s<ninvs; s++) {
      cout << "  " << symbols[s+1] << " = " << kin[0][s] << "," << endl;
    }
    printf("]\n");

  }
  fclose(card_fptr);

  //////
  // RUN DIRECTORY
  //////
  char *run_hash; 
  if (exit_sing == -1) {
    run_hash = strdup("amflow");
  } else {
    cout << endl; cout << "build run string..." << endl;
    char *run_str;
    point_string(
      &run_str,
      ninvs, symbols+1, kin[0]
    );
    cout << "run string: " << run_str << endl;
    calculate_hash(&run_hash, run_str);
    cout << "generate hash..." << endl;
  }
	cout << "run hash: " << run_hash << endl;

	if (dir_run) {
    join_path(&dir_run, dir_target, dir_run);
  } else {
    join_path(&dir_run, dir_target, run_hash);
	  join_path(&dir_run, dir_run, (char*)"/");
  }
  cout << endl; cout << "run directory: " << dir_run << endl;
  
  // create run directory
	if (make_dir(dir_run)) {
		printf("hash run directory already exists.\n");
	}

	//////
	// PREV RUN BOUNDARY
	//////
  int prev_run_bound = 0;
  char *dir_prev_run = NULL;
  if (prev_run_start_pt) {
    cout << endl; cout << "build prev-run string" << endl;
    char *prev_run_str;
    target_string(
      &prev_run_str,
      2, (char**)param_keys, param_values,
      ninvs, symbols+1, kin[0]
    );
    cout << "prev-run string: " << prev_run_str << endl;

    cout << "generate hash..." << endl;
    char *prev_run_hash;
    calculate_hash(&prev_run_hash, prev_run_str);
    cout << "target hash: " << prev_run_hash << endl;
    
    join_path(&dir_prev_run, dir_parent, (char*)"points/");
    join_path(&dir_prev_run, dir_prev_run, prev_run_hash);
    join_path(&dir_prev_run, dir_prev_run, (char*)"/");
    prev_run_bound = directory_exists(dir_prev_run);
    if (prev_run_bound) {
      cout << "prev run directory exists." << endl;
    }
  }

  // USE INPUT BOUNDARY
  int input_bound;
  card_fptr = fopen(filepath_card, "r");
  if (!read_param_int(&input_bound, (char*)"input-bound", card_fptr)) {
    input_bound = 0;  // default
  }
  fclose(card_fptr);

  //////
  // CACHE FILE PATHS
  //////
  char *dir_input = NULL, *dir_partial = NULL, *dir_output = NULL;
  join_path(&dir_input, dir_run, (char*)"input/");
  join_path(&dir_partial, dir_run, (char*)"partial/");
  join_path(&dir_output, dir_run, (char*)"output/");

  cout << endl; cout << "CACHE FOLDERS:" << endl;
  cout << dir_run << endl;
  cout << dir_input << endl;
  cout << dir_output << endl;
  cout << dir_partial << endl;

	make_dir(dir_input);
	make_dir(dir_partial);
	make_dir(dir_output);

  char *filepath_param = NULL, *filepath_PS = NULL, *filepath_bound = NULL;
  char *filepath_matrix = NULL, *filepath_roots = NULL, *filepath_branch_sing_lab = NULL;
  char *filepath_path = NULL, *filepath_path_tags = NULL, *filepath_sing_lab = NULL;
  char *filepath_log = NULL, *filepath_sol = NULL, *filepath_sol_interp = NULL;
  char *filepath_prev_run_PS = NULL; //, *filepath_prev_run_sol;
  char *filepath_target_info = NULL, *filepath_run_info = NULL;
  join_path(&filepath_target_info, dir_target, (char*)"info.txt");
  join_path(&filepath_run_info, dir_run, (char*)"info.txt");
  join_path(&filepath_param, dir_input, (char*)"parameters.txt");
  join_path(&filepath_PS, dir_input, (char*)"PS_points.txt");
  if (dir_bound) {
    if (dir_bound[strlen(dir_bound) -1] != '/') join_path(&dir_bound, dir_bound, (char*)"/");
    join_path(&filepath_bound, dir_bound, (char*)"bound");
  } else if (input_bound) {
    join_path(&filepath_bound, dir_common, (char*)"boundary/bound");
    // copy boundary into cache folder
    join_path(&tmp_filepath, dir_input, (char*)"boundary/");
    make_dir(tmp_filepath);
    join_path(&tmp_filepath, tmp_filepath, (char*)"bound");
    for (int ep=0; ep<eps_num; ep++) {
      snprintf(tmp_filepath1, MAX_PATH_LEN, "%s%d%s", filepath_bound, ep, file_ext);
      snprintf(tmp_filepath2, MAX_PATH_LEN, "%s%d%s", tmp_filepath, ep, file_ext);
      cout << "source file: " << tmp_filepath1 << endl;
      cout << "destin file: " << tmp_filepath2 << endl;
      copy_file(tmp_filepath2, tmp_filepath1);
    }
  } else {
    if (prev_run_bound) {
      // use previous run solutions as boundary
      join_path(&filepath_bound, dir_prev_run, (char*)"sol/sol");
    } else {
      join_path(&filepath_bound, dir_input, (char*)"boundary/bound");
    }
  }
  join_path(&filepath_matrix, dir_partial, (char*)"pfmat");
  join_path(&filepath_roots, dir_partial, (char*)"roots");
  join_path(&filepath_branch_sing_lab, dir_partial, (char*)"branch_sing_lab.txt");
  join_path(&filepath_path, dir_partial, (char*)"path");
  join_path(&filepath_path_tags, dir_partial, (char*)"path_tags");
  join_path(&filepath_sing_lab, dir_partial, (char*)"sing_lab");
  join_path(&filepath_log, dir_output, (char*)"log.txt");
  join_path(&filepath_sol, dir_output, (char*)"sol");
  join_path(&filepath_sol_interp, dir_output, (char*)"sol_interp.txt");

  cout << "prev run solutions: " << filepath_bound << endl;

  // SOLUTION PATH
  char *dir_outer_sol = NULL;
  char *filepath_outer_sol= NULL, *filepath_outer_sol_interp = NULL;
  join_path(&dir_outer_sol, dir_target, (char*)"sol/");
  cout << "outer sol directory: " << dir_outer_sol << endl;
  if (make_dir(dir_outer_sol)) {
    printf("hash outer sol directory already exists.\n");
  }
  join_path(&filepath_outer_sol, dir_outer_sol, (char*)"sol");
  join_path(&filepath_outer_sol_interp, dir_outer_sol, (char*)"sol_interp.txt");

  FILE* logfptr = stdout;

  // CACHE INFO
  if (opt_write > 0 || opt_write == -2) {
    log_params(2, (char**)param_keys, param_values, filepath_target_info, "w");
    log_point(ninvs, symbols+1, kin[1], (char*)"point", filepath_target_info, "a");
    
    if (exit_sing == -1) {
      join_path(&filepath_target_info, dir_amflow, (char*)"info.txt");
      log_params(2, (char**)param_keys, param_values, filepath_target_info, "w");
      log_point(ninvs, symbols+1, kin[1], (char*)"point", filepath_target_info, "a");
    } else {
      log_point(ninvs, symbols+1, kin[0], (char*)"point", filepath_run_info, "w");
    }
  }
  
  //////
  // BUILD PHASE-SPACE LINE
  //////
  if (exit_sing == -1) {
    // overwrite symbol array
    ninvs = 1;
    symbols[1] = strdup("eta");
    linearize_mass = 0;
    nmass = 0;
    is_mass[0] = 0;

    // overwrite phase-space strings
    kin[0][0] = strdup("0");
    kin[1][0] = strdup("1");
  }

	// take phase-space points as input
  mpc_t *PS_ini, *PS_fin;
  PS_ini = new mpc_t[ninvs];
  PS_fin = new mpc_t[ninvs];

  str_rk1_to_mpc_rk1(PS_ini, kin[0], ninvs);
  str_rk1_to_mpc_rk1(PS_fin, kin[1], ninvs);
  int *skip_inv = new int[ninvs];
  detect_active_invs(skip_inv, ninvs, PS_ini, PS_fin);

  // cout << endl; cout << "PHASE-SPACE LINE:" << endl;
  // cout << "initial point:" << endl;
  // for (int s=0; s<ninvs; s++) {
  //   print_mpc(&PS_ini[s]); cout << endl;
  // }
  // cout << "final point:" << endl;
  // for (int s=0; s<ninvs; s++) {
  //   print_mpc(&PS_fin[s]); cout << endl;
  // }
  // cout << endl; cout << "active invariants:" << endl;
  // for (int s=0; s<ninvs; s++) {
  //   cout << "s: " << s << ", skip: " << skip_inv[s] << endl;
  // }

  // LINEARIZE MASSES
  linearize_masses_mpc(
    PS_ini, PS_fin,
    nmass, mass
  );
  // cout << endl; cout << "linearized masses:" << endl;
  // cout << "initial point:" << endl;
  // for (int s=0; s<ninvs; s++) {
  //   print_mpc(&PS_ini[s]); cout << endl;
  // }
  // cout << "final point:" << endl;
  // for (int s=0; s<ninvs; s++) {
  //   print_mpc(&PS_fin[s]); cout << endl;
  // }

	//////
	// BUILD PHASE-SPACE POLY_FRAC
	//////

	// leave space for epsilon at the beginning of kin string
	char ***ep_kin;
	malloc_rk2_tens(ep_kin, 2, ninvs+1);
	for (int s=0; s<ninvs; s++) {
		ep_kin[0][s+1] = kin[0][s];
		ep_kin[1][s+1] = kin[1][s];
	}

	struct poly_frac *pspf = new struct poly_frac[ninvs+1];
	// cout <<endl; cout << "PS pf:" << endl;
	for (int s=1; s<=ninvs; s++) {
		// cout << "s = " << s << endl;
		// cout << "skip_inv, is_mass = " << skip_inv[s-1] << ", " << is_mass[s-1] << endl;
		poly_frac_build(&pspf[s]);
		generate_PS_pf(&pspf[s], ep_kin, skip_inv[s-1], is_mass[s-1], s, 1, wp2);
		// poly_frac_print(&pspf[s]);
	}

  //////
  // INPUT MATRICES
  //////
  // INPUT MATRICES AS STRINGS
	char ****mats_str = NULL;
  if (exit_sing != -1) {
  // if (1) {
    cout << "load DE matrices from " << filepath_mats << endl;
    cout << endl; cout << "INPUT MATRICES:" << endl;
    mats_str = new char***[ninvs];
    load_DE_matrices(
      mats_str, &dim,
      ninvs, symbols, skip_inv,
      filepath_mats, file_ext
    );

  }

  cout << "dim = " << dim << endl;

  int *starting_ord = new int[dim];

  //////
  // PREPARATIONS
  //////
  // ROOTS
  int *zero_label = new int[eps_num];
  int *nroots = new int[eps_num];
  mpc_t **roots = new mpc_t*[eps_num];
  mpfr_t **tols = new mpfr_t*[eps_num];
  int **perm = new int*[eps_num];

  // DE
  poly_frac ***pfmat;
  malloc_rk3_tens(pfmat, eps_num, dim, dim);
  poly_frac_rk3_build(pfmat, eps_num, dim, dim);

  // PATH
  mpc_t **path = new mpc_t*[eps_num];
  int **path_tags = new int*[eps_num];
  int *neta_values = new int[eps_num];
  int  **sing_lab = new int*[eps_num];
  int *nsings = new int[eps_num];
  mpc_t **path_mp = new mpc_t*[eps_num];
  int **path_tags_mp = new int*[eps_num];
  int *neta_values_mp = new int[eps_num];
  int  **sing_lab_mp = new int*[eps_num];
  int *nsings_mp = new int[eps_num];

  // SOLUTIONS
  mpc_t ***solutions;
  solutions = new mpc_t**[1];
  malloc_rk2_tens(*solutions, dim, eta_ord+1);
  init_rk2_mpc(*solutions, dim, eta_ord+1);
  mpc_t **sol_at_eps;
  malloc_rk2_tens(sol_at_eps, dim, eps_num);
  init_rk2_mpc(sol_at_eps, dim, eps_num);

  // WRT/CMP
  mpc_t **sol_at_eps_wrt_cmp;
  int dim_wrt_cmp;

  //////
  // GENERATE BOUNDARIES
  //////
  card_fptr = fopen(filepath_card, "r");
  int gen_bound;
  if (!read_param_int(&gen_bound, (char*)"gen-bound", card_fptr)) {
    cout << "gen-bound option not specified" << endl;
    gen_bound = 0;
  }
  cout << endl; cout << "gen-bound: " << gen_bound << endl;
  if (gen_bound == 1 && exit_sing == 0) {
    perror("boundary generation activated with exit-sing = 0");
    exit(1);
  }
  if (gen_bound == 1 && input_bound) {
    perror("boundary generation activated but bound-dir also specified");
    exit(1);
  }

  mpc_t **bound = NULL; // NULL needed
  if (gen_bound == 1) {
    cout << endl; cout << "GENERATE BOUNDARIES" << endl;
    malloc_rk2_tens(bound, eps_num, dim);
    init_rk2_mpc(bound, eps_num, dim);
    generate_boundaries(
      bound,
      filepath_bound_build, eps_num, dim, nloops,
      eps_str, kin, symbols, ninvs
    );

    cout << "generated boundaries (1st epsilon):" << endl;
    for (int ep=0, i=0; i<dim; i++) {
      cout << "BC n." << i << ": "; print_mpc(&bound[ep][i]); cout << endl;
    }
  }

  //////
  // ANALYSE BRANCH CUTS
  //////
	int nbranch_roots = 1;
	mpc_t *branch_roots = new mpc_t[1];
	mpc_init3(branch_roots[0], wp2, wp2);
	mpc_set_ui(branch_roots[0], 0, MPFR_RNDN);
	mpfr_t *branch_tols = new mpfr_t[1];
	mpfr_init2(branch_tols[0], wp2);
	mpfr_set(branch_tols[0], mpfr_tol, MPFR_RNDN);

  int num_branch_roots_tmp;
  int **branch_sing_lab;
  malloc_rk2_tens(branch_sing_lab, nbranches, 3);
  int *branch_deg = new int[nbranches];
  mpc_t **branch_poly = new mpc_t*[nbranches];
  mpc_t *branch_roots_tmp = new mpc_t[2*nbranches];
  mpfr_t *branch_tols_tmp = new mpfr_t[2*nbranches];

  process_branch_points(
    &num_branch_roots_tmp, branch_roots_tmp, branch_tols_tmp,
    branch_sing_lab,
    branch_deg, branch_poly,
    &nbranch_roots, &branch_roots, &branch_tols,
    nbranches, branch_cut_str,
    ninvs, symbols, is_mass,
    skip_inv, ep_kin, pspf
  );
  // cout << endl; cout << num_branch_roots_tmp << " branch roots:" << endl;
  // print_poly(branch_roots_tmp, num_branch_roots_tmp-1);
  // cout << "branch sing labels:" << endl;
  // for (int b=0; b<nbranches; b++) {
  //   for (int k=1; k<=branch_deg[b]; k++) {
  //     cout << branch_sing_lab[b][k] << endl;
  //   }
  // }
  // cout << endl; cout << nbranch_roots << " roots:" << endl;
  // print_poly(branch_roots, nbranch_roots-1);

  //////
  // LOAD BOUNDARY
  //////

  //////
  // GO TO CHECKPOINT
  //////
  switch (opt_checkpoint) {
    case -2:
    case -1:
    case 0:
      break;
    case 1:
      goto goto_checkpoint1;
    case 2:
      goto goto_checkpoint2;
    default:
      perror("no valid checkpoint option");
      exit(1);
  }

  //////
  // GENERATE DE
  //////
  if (exit_sing == -1) {
  // if (0 && exit_sing == -1) {
    kira_to_DE_pf(
      pfmat,
      zero_label, nroots, roots, tols,
      nbranch_roots, branch_roots, branch_tols,
      ninvs, symbols, is_mass,
      pspf,
      skip_inv, ep_kin,
      dim, mats_str,
      nbranches, branch_deg,
      eps_num, eps_str,
      MI_eta, dir_amflow,
      terminal
    );

  } else {
    generate_poly_frac_DE(
      pfmat,
      zero_label, nroots, roots, tols,
      nbranch_roots, branch_roots, branch_tols,
      ninvs, symbols, is_mass,
      pspf,
      skip_inv, ep_kin,
      dim, mats_str,
      nbranches, branch_deg,
      eps_num, eps_str,
      terminal
    );
  }

  if (print) {
  cout << endl; cout << "DE MATRIX (1st epsilon):" << endl;
  poly_frac_rk2_print(pfmat[0], dim, dim);
  }
  
  cout << endl; cout << "ROOTS (1st epsilon):" << endl;
  cout << "zero_label = " << zero_label[0] << endl;
  cout << nroots[0] << " roots" << endl;
  print_poly(roots[0], nroots[0]-1);

  //////
  // BUILD PATH
  //////
  cout << endl; cout << "BUILD PATH..." << endl;
  for (int ep=0; ep<eps_num; ep++) {
    if (exit_sing == -1) {
      get_path_PS_infty_mp(
        &path[ep], &path_tags[ep], &neta_values[ep], &sing_lab[ep], &nsings[ep],
        roots[ep], nroots[ep], zero_label[ep]
      );
    } else {
      get_path_PS_mp(
        &path[ep], &path_tags[ep], &neta_values[ep], &sing_lab[ep], &nsings[ep],
        roots[ep], nroots[ep], zero_label[ep]
      );
    }
  }

  cout << endl; cout << "SINGULAR POINTS (1st epsilon):" << endl;
  cout << nsings[0] << " points" << endl;
  for (int i=0; i<nsings[0]; i++) {
    cout << "lab: " << sing_lab[0][i] << ", root: "; print_mpc(&roots[0][sing_lab[0][i]]); cout << endl;
  }
  cout << endl; cout << "PATH (1st epsilon):" << endl;
  cout << neta_values[0] << " points" << endl;
  for (int i=0; i<neta_values[0]; i++) {
    cout << i << ". tag: " << path_tags[0][i] << ", point: ";
    print_mpc(&path[0][i]); cout << endl;
  }

  //////
  // CHECKPOINT 1
  //////
  int write;
  switch (opt_write) {
    case -3:
    case -1:
    case 2:
      write = -1;
      break;
    case 0:
    case 1:
      write = opt_write;
      break;
    case -2:
      write = 1;
      break;
  }

  if (write >= 0) {
    wrt_cmp_DE(
      perm,
      zero_label, nroots, roots,
      pfmat, eps_num, dim,
      0.80,
      file_ext, filepath_matrix, filepath_roots, opt_write
    );
    wrt_cmp_path(
      path, path_tags, eps_num, neta_values,
      nsings, sing_lab, perm,
      file_ext, filepath_path, filepath_path_tags, filepath_sing_lab, opt_write
    );
  }

  goto_checkpoint1:
  switch (opt_checkpoint) {
    case -1:
      goto goto_exit;
    case 1:
      //////
      // LOAD POLY_FRAC DE FROM FILE
      //////
      for (int ep=0; ep<eps_num; ep++) {
        // MATRIX
        snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_matrix, ep, file_ext);
        cout << endl; cout << "reading poly_frac DE from file " << tmp_filepath << endl;
        poly_frac_rk2_from_file(tmp_filepath, pfmat[ep], dim, dim);

        // ROOTS
        snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_roots, ep, file_ext);
        cout << "reading mpc_t roots from file " << tmp_filepath << endl;
        nroots[ep] = count_lines(tmp_filepath) - 1;
        // if (roots[ep]) {
        //   delete[] roots[ep];
        // }
        roots[ep] = new mpc_t[nroots[ep]];
        init_rk1_mpc(roots[ep], nroots[ep]);
        int_rk0_mpc_rk1_from_file(tmp_filepath, roots[ep], nroots[ep], &zero_label[ep]);
        // snprintf(tmp_filepath, sizeof(tmp_filepath), "%s%d%s", filepath_branch_sing_lab, ep, file_ext);
        // cout << endl; cout << "read branch sing labs from file " << tmp_filepath << endl;
        // int_rk2_from_file(tmp_filepath, branch_sing_lab, 3, 2);

        // PATH
        // points
        snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_path, ep, file_ext);
        cout << "reading path points from file " << tmp_filepath << endl;
        neta_values[ep] = count_lines(tmp_filepath);
        path[ep] = new mpc_t[neta_values[ep]];
        init_rk1_mpc(path[ep], neta_values[ep]);
        mpc_rk1_from_file(tmp_filepath, path[ep], neta_values[ep]);
        // tags
        snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_path_tags, ep, file_ext);
        cout << "reading path tags from file " << tmp_filepath << endl;
        path_tags[ep] = new int[neta_values[ep]];
        int_rk1_from_file(tmp_filepath, path_tags[ep], neta_values[ep]);
        // singular labels
        snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_sing_lab, ep, file_ext);
        cout << "reading singular labels from file " << tmp_filepath << endl;
        nsings[ep] = count_lines(tmp_filepath);
        sing_lab[ep] = new int[nsings[ep]];
        int_rk1_from_file(tmp_filepath, sing_lab[ep], nsings[ep]);
      }

      cout << endl; cout << "ROOTS (1st epsilon):" << endl;
      cout << "zero_label = " << zero_label[0] << endl;
      cout << nroots[0] << " roots:" << endl;
      print_poly(roots[0], nroots[0]-1);
      // cout << endl; cout << "ROOTS (2nd epsilon):" << endl;
      // cout << "zero_label = " << zero_label[1] << endl;
      // cout << nroots[1] << " roots:" << endl;
      // print_poly(roots[1], nroots[1]-1);

      cout << endl; cout << "SINGULAR POINTS (1st epsilon):" << endl;
      cout << nsings[0] << " points" << endl;
      for (int i=0; i<nsings[0]; i++) {
        cout << "lab: " << sing_lab[0][i] << ", root: "; print_mpc(&roots[0][sing_lab[0][i]]); cout << endl;
      }
      cout << endl; cout << "PATH (1st epsilon):" << endl;
      cout << neta_values[0] << " points" << endl;
      for (int i=0; i<neta_values[0]; i++) {
        cout << i << ". tag: " << path_tags[0][i] << ", point: ";
        print_mpc(&path[0][i]); cout << endl;
      }

      // // MP PATH
      // cout << endl; cout << "BUILD PATH..." << endl;
      // for (int ep=0; ep<eps_num; ep++) {
      //   if (exit_sing == -1) {
      //     get_path_PS_infty_mp(
      //       &path_mp[ep], &path_tags_mp[ep], &neta_values_mp[ep], &sing_lab_mp[ep], &nsings_mp[ep],
      //       roots[ep], nroots[ep], zero_label[ep]
      //     );
      //   } else {
      //     get_path_PS_mp(
      //       &path_mp[ep], &path_tags_mp[ep], &neta_values_mp[ep], &sing_lab_mp[ep], &nsings_mp[ep],
      //       roots[ep], nroots[ep], zero_label[ep]
      //     );
      //   }
      // }
      // MP PATH
      // cout << endl; cout << "BUILD PATH..." << endl;
      // for (int ep=0; ep<eps_num; ep++) {
      //   if (exit_sing == -1) {
      //     get_path_PS_infty_mp(
      //       &path[ep], &path_tags[ep], &neta_values[ep], &sing_lab[ep], &nsings[ep],
      //       roots[ep], nroots[ep], zero_label[ep]
      //     );
      //   } else {
      //     get_path_PS(
      //       &path[ep], &path_tags[ep], &neta_values[ep], &sing_lab[ep], &nsings[ep],
      //       roots[ep], nroots[ep], zero_label[ep]
      //     );
      //   }
      // }

      // cout << endl; cout << "SINGULAR POINTS MP (1st epsilon):" << endl;
      // cout << nsings_mp[0] << " points" << endl;
      // for (int i=0; i<nsings_mp[0]; i++) {
      //   cout << "lab: " << sing_lab_mp[0][i] << ", root: "; print_mpc(&roots[0][sing_lab_mp[0][i]]); cout << endl;
      // }
      // cout << endl; cout << "PATH MP (1st epsilon):" << endl;
      // cout << neta_values_mp[0] << " points" << endl;
      // for (int i=0; i<neta_values_mp[0]; i++) {
      //   cout << i << ". tag: " << path_tags_mp[0][i] << ", point: ";
      //   print_mpc(&path_mp[0][i]); cout << endl;
      // }

      // wrt_cmp_DE(
      //   perm,
      //   zero_label, nroots, roots,
      //   pfmat, eps_num, dim,
      //   0.80,
      //   file_ext, filepath_matrix, filepath_roots, 0
      // );
      // wrt_cmp_path(
      //   path_mp, path_tags_mp, eps_num, neta_values_mp,
      //   nsings_mp, sing_lab_mp, perm,
      //   file_ext, filepath_path, filepath_path_tags, filepath_sing_lab, 1
      // );
      break;
  
    case -2:
    case 0:
      break;
    default:
      perror("no valid checkpoint option");
      exit(1);
  }

  // cache path info
  if (opt_write > 0 || opt_write == -2) {
    log_param((char*)"num-propag", neta_values[0]-1, filepath_run_info, "a");
    log_param((char*)"num-sing-propag", nsings[0], filepath_run_info, "a");
  }

  //////
  // SOLVE DE
  //////

  //////
  // EPSILON LOOP
  //////
  // goto goto_interpolate;
  if (opt_all_eps == 0) {
    eps_num = 1;
  }
  propagate_all_eps(
    // OUTPUT
    sol_at_eps,
    // INPUT
    exit_sing, nloops,
    pfmat,
    zero_label, nroots, roots,
    neta_values, path, path_tags, nsings, sing_lab,
    solutions, dim, eta_ord,
    eps_num, eps_str,
    ninvs, PS_ini, PS_fin, symbols,
    is_mass, skip_inv,
    nbranches, branch_deg, branch_poly, branch_sing_lab,
    // bound_behav, mi_eig, mi_eig_num,
    gen_bound, filepath_bound, bound,
    filepath_bound_build, filepath_bound_behav,
    // filepath_matrix, filepath_roots, // filepath_branch_sing_lab,
    filepath_path, filepath_path_tags, filepath_sol,
    file_ext, logfptr, opt_write,
    terminal
  );

  if (exit_sing == -1) {
    malloc_rk2_tens(sol_at_eps_wrt_cmp, dim_eta_less, eps_num);
    init_rk2_mpc(sol_at_eps_wrt_cmp, dim_eta_less, eps_num);
    //////
    // SELECT ETA-LESS MIs FROM SOLUTION
    //////
    cout << endl; cout << "SELECT ETA-LESS MIs FROM SOLUTION" << endl;
    for (int m=0; m<dim_eta_less; m++) {
      // cout << "m, MI_idx = " << m << ", " << MI_idx[m] << endl;
      for (int ep=0; ep<eps_num; ep++) {
        mpc_set(sol_at_eps_wrt_cmp[m][ep], sol_at_eps[MI_idx[m]][ep], MPFR_RNDN);
        // print_mpc(&sol_at_eps_wrt_cmp[m][ep]); cout << endl;
      }
    }
    dim_wrt_cmp = dim_eta_less;
  } else {
    sol_at_eps_wrt_cmp = sol_at_eps;
    dim_wrt_cmp = dim;
  }

  //////
  // WRITE/COMPARE RESULTS
  //////
  if (opt_write > 0 || opt_write == -2) {
    // cache solution
    join_path(&tmp_filepath, filepath_sol, file_ext);
    mpc_rk2_to_file(tmp_filepath, sol_at_eps, dim, eps_num);
  }

  if (opt_write > 0)  {
    // outer solution
    // #pair
    // per-epsilon solutions
    for (int ep=0; ep<eps_num; ep++) {
      snprintf(tmp_filepath, MAX_PATH_LEN, "%s%d%s", filepath_outer_sol, ep, file_ext);
      FILE *sol_fptr = fopen(tmp_filepath, "w");
      for (int i=0; i<dim_wrt_cmp; i++) {
        mpc_out_str(sol_fptr, 10, 0, sol_at_eps_wrt_cmp[i][ep], MPFR_RNDN); fprintf(sol_fptr, "\n");
      }
      fclose(sol_fptr);
    }
    // all-epsilon solutions
    join_path(&tmp_filepath, filepath_outer_sol, file_ext);
    mpc_rk2_to_file(tmp_filepath, sol_at_eps_wrt_cmp, dim_wrt_cmp, eps_num);

  } else if (opt_write != -3) {
    // enlarge tolerance
    mpfr_t mpfr_tol_orig;
    mpfr_init2(mpfr_tol_orig, wp2);
    mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
    mpfr_tol_set_wp(precision*2);
    // cout << endl; cout << "check with enlarged mpfr tol:" << endl;
    // mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;
    
    join_path(&tmp_filepath, filepath_outer_sol, file_ext);
    cout << endl; cout << "reading from " << tmp_filepath << endl;
    int bench_eps_num = count_lines(tmp_filepath)/dim_wrt_cmp;
    mpc_t **bench_sol_at_eps;
    malloc_rk2_tens(bench_sol_at_eps, dim_wrt_cmp, bench_eps_num);
    init_rk2_mpc(bench_sol_at_eps, dim_wrt_cmp, bench_eps_num);
    mpc_rk2_from_file(tmp_filepath, bench_sol_at_eps, dim_wrt_cmp, bench_eps_num);
    cout << "CHECK result..." << endl;
    mpc_t **bench_sol_at_eps_sliced;
    malloc_rk2_tens(bench_sol_at_eps_sliced, dim_wrt_cmp, eps_num);
    init_rk2_mpc(bench_sol_at_eps_sliced, dim_wrt_cmp, eps_num);
    for (int m=0; m<dim_wrt_cmp; m++) {
      for (int ep=0; ep<eps_num; ep++) {
        mpc_set(bench_sol_at_eps_sliced[m][ep], bench_sol_at_eps[m][ep], MPFR_RNDN);
      }
    }
    mpc_rk2_compare(bench_sol_at_eps, sol_at_eps_wrt_cmp, dim_wrt_cmp, eps_num);

    // restore original tol
    mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
  }

  if (!opt_all_eps) {
    cout << "processing only one epsilon: skip interpolatation" << endl;
    goto goto_exit;
  }

  goto_checkpoint2:
  switch(opt_checkpoint) {
    case 0:
    case 1:
      break;
    case -2:
      goto goto_exit;
    case 2:
      // LOAD FROM FILE
      join_path(&tmp_filepath, filepath_sol, file_ext);
      mpc_rk2_from_file(tmp_filepath, sol_at_eps, dim, eps_num);
      cout << "loaded results (# MI) X (# EPS):" << endl;
      print_rk2_mpc(sol_at_eps, dim, eps_num);
      break;
    default:
      perror("no valid checkpoint option");
      exit(1);
  }

  //////
  // INTERPOLATE EPSILON ORDERS
  //////
  cout << endl; cout << "INTERPOLATE EPSILON ORDERS..." << endl;
  mpc_t **sol_eps_ord;
  malloc_rk2_tens(sol_eps_ord, dim, eps_num);
  init_rk2_mpc(sol_eps_ord, dim, eps_num);
  mpc_t **sol_eps_ord_wrt_cmp;
  interpolate_epsilon_orders(
    sol_eps_ord,
    sol_at_eps, eps_str,
    dim, eps_num, nloops,
    precision,
    NULL
  );

  fprintf(logfptr, "\nEPSILON ORDERS:\n");
  for (int i=0; i<dim; i++) {
    fprintf(logfptr, "MI n. %d\n", i);
    for (int ep=0; ep<=order+2*nloops; ep++) {
      fprintf(logfptr, "eps^%d: ", ep - 2*nloops);
      // fprintf(logfptr, "eps order %d\n", ep - 2*nloops + starting_ord[i]);
      mpc_out_str(logfptr, 10, 0, sol_eps_ord[i][ep], MPFR_RNDN); fprintf(logfptr, "\n");
    }
  }

  mpc_t **sol_eps_ord_pruned;
  if (opt_prune_eps_mode) {
    malloc_rk2_tens(sol_eps_ord_pruned, dim, eps_num);
    init_rk2_mpc(sol_eps_ord_pruned, dim, eps_num);
    if (opt_prune_eps_mode == 1) {
      // RELATIVE PRUNE
      cout << endl; cout << "RELATIVE EPSILON PRUNE" << endl;
      interpolate_epsilon_orders_prune(
        sol_eps_ord_pruned,
        sol_at_eps, eps_str,
        dim, eps_num, nloops,
        precision, order,
        starting_ord
      );

      for (int m=0; m<dim; m++) {
        if (starting_ord[m] == 0) {
          continue;
        }

        for (int ep=eps_num-1; ep>=starting_ord[m]; ep--) {
          mpc_set(sol_eps_ord_pruned[m][ep], sol_eps_ord_pruned[m][ep-starting_ord[m]], MPFR_RNDN);
        }
        for (int ep=0; ep<starting_ord[m]; ep++) {
          mpc_set_ui(sol_eps_ord_pruned[m][ep], 0, MPFR_RNDN);
        }
      }

    } else {
      // ABSOLUTE PRUNE
      cout << endl; cout << "ABSOLUTE EPSILON PRUNE" << endl;
      // wp2 = (precision+1) * log2(10);
      for (int i=0; i<dim; i++) {
        for (int ep=0; ep<=order+2*nloops; ep++) {
          // prune real part
          if (mpfr_get_exp(mpc_realref(sol_eps_ord[i][ep])) < prune_eps_exp2) {
            mpfr_set_ui(mpc_realref(sol_eps_ord_pruned[i][ep]), 0, MPFR_RNDN);
          } else {
            mpfr_set(mpc_realref(sol_eps_ord_pruned[i][ep]), mpc_realref(sol_eps_ord[i][ep]), MPFR_RNDN);
          }
          // prune imag part
          if (mpfr_get_exp(mpc_imagref(sol_eps_ord[i][ep])) < prune_eps_exp2) {
            mpfr_set_ui(mpc_imagref(sol_eps_ord_pruned[i][ep]), 0, MPFR_RNDN);
          } else {
            mpfr_set(mpc_imagref(sol_eps_ord_pruned[i][ep]), mpc_imagref(sol_eps_ord[i][ep]), MPFR_RNDN);
          }
        }
      }
    }
  } else {
    // NO PRUNE
    sol_eps_ord_pruned = sol_eps_ord;
  }

  // get minimum exponent
  int min_exp;
  min_exp = mpfr_get_exp(mpc_realref(sol_eps_ord_pruned[0][0]));
  for (int m=0; m<dim; m++) {
    for (int ep=0; ep<order+2*nloops; ep++) {
      if (mpfr_get_exp(mpc_realref(sol_eps_ord_pruned[m][ep])) < min_exp) {
        min_exp = mpfr_get_exp(mpc_realref(sol_eps_ord_pruned[m][ep]));
      }
      if (mpfr_get_exp(mpc_imagref(sol_eps_ord_pruned[m][ep])) < min_exp) {
        min_exp = mpfr_get_exp(mpc_imagref(sol_eps_ord_pruned[m][ep]));
      }
    }
  }
  min_exp = log10(2) * min_exp;
  // fprintf(stderr, "MINIMUM EXPONENT: %d\n", min_exp);
  // fprintf(stdout, "MINIMUM EXPONENT: %d\n", min_exp);

  if (exit_sing == -1) {
    malloc_rk2_tens(sol_eps_ord_wrt_cmp, dim_eta_less, eps_num);
    init_rk2_mpc(sol_eps_ord_wrt_cmp, dim_eta_less, eps_num);
    //////
    // SELECT ETA-LESS MIs FROM SOLUTION
    //////
    for (int m=0; m<dim_eta_less; m++) {
      for (int ep=0; ep<eps_num; ep++) {
        // mpc_set(sol_eps_ord_wrt_cmp[m][ep], sol_eps_ord[MI_idx[m]][ep], MPFR_RNDN);
        mpc_set(sol_eps_ord_wrt_cmp[m][ep], sol_eps_ord_pruned[MI_idx[m]][ep], MPFR_RNDN);
      }
    }
    dim_wrt_cmp = dim_eta_less;

  } else {
    sol_eps_ord_wrt_cmp = sol_eps_ord_pruned;
    dim_wrt_cmp = dim;
  }

  // if (opt_prune_eps_mode || exit_sing == -1) {
  //   fprintf(logfptr, "\nEPSILON ORDERS (");
  //   if (opt_prune_eps_mode) fprintf(logfptr, "pruned, ");
  //   if (exit_sing == -1) fprintf(logfptr, "eta-less, ");
  //   fprintf(logfptr, "):\n");
  //   for (int i=0; i<dim_wrt_cmp; i++) {
  //     fprintf(logfptr, "MI n. %d\n", i);
  //     for (int ep=0; ep<=order+2*nloops; ep++) {
  //       fprintf(logfptr, "eps^%d: ", ep - 2*nloops);
  //       // fprintf(logfptr, "eps order %d\n", ep - 2*nloops + starting_ord[i]);
  //       mpc_out_str(logfptr, 10, 0, sol_eps_ord_wrt_cmp[i][ep], MPFR_RNDN); fprintf(logfptr, "\n");
  //     }
  //   }
  // }

  if (opt_write > 0 || opt_write == -2) {
    // cache file
    cout << endl; cout << "writing to " << filepath_sol_interp << endl;
    mpc_rk2_to_file(filepath_sol_interp, sol_eps_ord, dim, order+2*nloops+1);
  }

  if (opt_write > 0) {
    // outer solutions
    cout << endl; cout << "writing to " << filepath_outer_sol_interp << endl;
    mpc_rk2_to_file(filepath_outer_sol_interp, sol_eps_ord_wrt_cmp, dim_wrt_cmp, order+2*nloops+1);
  } else if (opt_write != -3) {
    mpfr_t mpfr_tol_orig;
    mpfr_init2(mpfr_tol_orig, wp2);
    mpfr_set(mpfr_tol_orig, mpfr_tol, MPFR_RNDN);
    // mpfr_tol_set_wp(precision*2);
    mpfr_tol_set_wp(precision*1.2);
    // cout << endl; cout << "check with enlarged mpfr tol:" << endl;
    // mpfr_out_str(stdout, 10, 0, mpfr_tol, MPFR_RNDN); cout << endl;

    cout << endl; cout << "CHECK interpolation..." << endl;
    mpc_t **bench_sol_eps_ord;
    malloc_rk2_tens(bench_sol_eps_ord, dim_wrt_cmp, order+2*nloops+1);
    init_rk2_mpc(bench_sol_eps_ord, dim_wrt_cmp, order+2*nloops+1);
    // cout << "benchmark:" << endl;
    // print_rk2_mpc(bench_sol_eps_ord, dim_wrt_cmp, order+2*nloops+1);
    // cout << "current:" << endl;
    // print_rk2_mpc(sol_eps_ord_wrt_cmp, dim_wrt_cmp, order+2*nloops+1);
    mpc_rk2_from_file(filepath_outer_sol_interp, bench_sol_eps_ord, dim_wrt_cmp, order+2*nloops+1);
    // mpc_rk2_compare(bench_sol_eps_ord, sol_eps_ord_wrt_cmp, dim_wrt_cmp, order+2*nloops+1);
    mpc_rk2_compare(bench_sol_eps_ord, sol_eps_ord_wrt_cmp, dim_wrt_cmp, order+2*nloops+1);

    // restore original tol
    mpfr_set(mpfr_tol, mpfr_tol_orig, MPFR_RNDN);
  }

  // FINAL PRINT
  fprintf(logfptr, "\nEPSILON ORDERS (target precision):\n");
  print_result(logfptr, precision, sol_eps_ord_wrt_cmp, dim_wrt_cmp, order, nloops);

  if (filepath_result) {
    FILE *resfptr = fopen(filepath_result, "w");
    if (!resfptr) {
		  fprintf(stderr, "error while creating/opening file %s\n", filepath_result);
		  exit(1);
	  }
    // print_result(resfptr, precision, sol_eps_ord_prec, dim_wrt_cmp, order, nloops);
    print_result(resfptr, precision, sol_eps_ord_wrt_cmp, dim_wrt_cmp, order, nloops);
  }

  goto_exit:
  //////
  // EXIT PROGRAM
  //////
  // FREE
  delete[] zero_label;
  for (int ep=0; ep<eps_num; ep++) {
    mpc_rk1_clear(roots[ep], nroots[ep]);
    delete[] roots[ep];
    mpfr_rk1_clear(tols[ep], nroots[ep]);
    delete[] tols[ep];
  }
  delete[] roots;
  delete[] tols;
  delete[] nroots;
  poly_frac_rk3_free(pfmat, eps_num, dim, dim);
  del_rk3_tens(pfmat, eps_num, dim);
  mpc_rk2_clear(*solutions, dim, eta_ord+1);
  del_rk2_tens(*solutions, dim);  
  delete[] solutions;
  mpc_rk2_clear(sol_at_eps, dim, eps_num);
  del_rk2_tens(sol_at_eps, dim);

  // close progress bars
  fprintf(terminal, "\033[2K\r"); fflush(terminal); usleep(sleep_time);
  if (!opt_bar) {
    fclose(terminal);
  }
  return 0;

}
