/********************************
* MCCE library header file
********************************/

#include <stdio.h>

/*--- Constants ---*/
/* Constants for geometry transformation recorder */
#define USERERR -1
#define MAXCHAR_LINE 160
#define DUMMY_GDBM   "~temp.dbm.XXXXXX"
#define FN_RUNPRM    "run.prm"
#define FN_HVROT     "hvrot.pdb"
#define STEP1_OUT    "step1_out.pdb"
#define STEP2_OUT    "step2_out.pdb"
#define STEP3_OUT    "step3_out.pdb"
#define ENERGY_TABLE "energies.opp"
#define FN_CONFLIST1 "head1.lst"
#define FN_CONFLIST2 "head2.lst"
#define FN_CONFLIST3 "head3.lst"
#define ROTSTAT      "rot_stat"
#define MC_OUT       "mc_out"
#define DETAIL       "fort.36"
#define OCC_TABLE    "fort.38"
#define TOT_CRG      "sum_crg.out"
#define CURVE_FITTING "pK.out"
#define MFE_OUT       "analysis"
#define LIST_ROT_GOLD "list_rot.gold"
#define ROT_GOLD      12

/* at 298.15K */
#define ROOMT        298.15
#define KCAL2KT      1.688
#define PH2KCAL      1.364

#define  MAXNBORS     10
#define  N_ELEM_MAX     60
#define  VDW_CUTOFF_NEAR  1
#define  VDW_ELIMIT_NEAR  999
#define  VDW_CUTOFF_FAR   10
#define  VDW_ELIMIT_FAR   0

/*--- Global variables ---*/
typedef struct {
    char inpdb[256];
    char param[256];
    char mcce_home[256];
    char debug_log[256];
    char progress_log[256];
    char new_tpl[256];
    char extra[256];

    char do_premcce;
    char do_rotamers;
    char do_energies;
    char do_monte;

    char quick_energies;

    float default_radius;
    float factor_14lj;
    float epsilon_columb;

    char  terminals;
    float clash_distance;
    float h2o_sascutoff;
    char  rename_rules[256];

    char  rot_specif;
    char  rot_swap;
    char  swing;
    char  pack;

    float prune_thr;

    char  hdirected;
    float hdirdiff;
    int   hdirlimt;

    float sas_cutoff;
    float vdw_cutoff;
    float repack_cutoff;
    float ngh_vdw_thr;
    float repack_e_thr_exposed;
    float repack_e_thr_buried;
    float phi_swing;
    int   rotations;
    int   rotamer_limit;
    int   repacks;
    int   nconf_limit;

    char  relax_wat;
    float water_relax_thr;

    int   hv_relax_ncycle;
    int   hv_relax_niter;
    float hv_relax_vdw_thr;
    float hv_relax_hv_vdw_thr;
    float hv_relax_dt;
    float hv_tors_scale;
    float hv_relax_constraint;
    float hv_relax_constraint_frc;
    int   hv_relax_n_shake;
    float hv_relax_shake_tol;
    int   hv_relax_include_ngh;
    float hv_relax_ngh_thr;

    char  relax_h;
    float relax_e_thr;
    int   relax_nstates;
    int   relax_n_hyd;
    float relax_clash_thr;
    float relax_phi;
    int   relax_niter;
    float relax_torq_thr;
    float prune_rmsd;
    float prune_ele;
    float prune_vdw;

    char  reassign;
    float epsilon_prot;
    float epsilon_solv;
    int   grids_delphi;
    float grids_per_ang;
    float radius_probe;
    float ionrad;
    float salt;
    char  delphi_exe[256];
    char  delphi_fails;
    int   err_msg;
    char  recalc_tors;

    char  average_pairwise;
    float warn_pairwise;

    int   monte_seed;
    float monte_temp;
    int   monte_flips;
    int   monte_nstart;
    int   monte_neq;
    float monte_reduce;

    int   monte_runs;
    int   monte_niter;
    int   monte_trace;
    int   nstate_max;

    char  monte_adv_opt;
    char  adding_conf;
    char  monte_old_input;
    int   monte_n_red;
    int   monte_niter_cycle;
    int   monte_niter_min;
    int   monte_niter_max;
    int   monte_niter_chk;
    float monte_converge;
    int   monte_do_energy;
    int   monte_print_nonzero;

    float anneal_temp_start;
    float anneal_nstep;
    float anneal_niter_step;

    char  titr_type;
    float titr_ph0;
    float titr_phd;
    float titr_eh0;
    float titr_ehd;
    int   titr_steps;

    float big_pairwise;

    /* advanced features */
    char rotate_res[256];
    float vdwf1;
    float vdwf2;
    int   delphi_start;
    int   delphi_end;

    float scale_vdw0;
    float scale_vdw1;
    float scale_vdw;
    float scale_tor;
    float scale_ele;
    float scale_dsolv;

    int   skip_ele;
    char  delphi_folder[256];
    char  delphi_clean;
    float PI;
    float d2r;
} ENV;

extern ENV env;

/*--- Data prototype ---*/
/* data types for geometry objects */
typedef struct {
    float M[4][4];
} GEOM;

typedef struct {
    float x;
    float y;
    float z;
} VECTOR;

typedef struct {
    VECTOR p0;                 /* line through this point */
    VECTOR t;                  /* direction cosine of the line */
} LINE;

typedef struct {
    VECTOR p0;                 /* plane through this point */
    VECTOR t;                  /* direction cosine of the plane */
} PLANE;


/*--- String structure types ---*/
typedef struct {              /* an array of n strings */
    int  n;
    char **strings;
} STRINGS;


/*--- Geometry functions ---*/
void mxm4(float m1[4][4], float m2[4][4], float m3[4][4]);
void geom_reset(GEOM *op);
void geom_move(GEOM *op, VECTOR v);
void geom_roll(GEOM *op, float phi, LINE axis);
GEOM geom_3v_onto_3v(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR t1, VECTOR t2, VECTOR t3);
void geom_inverse(GEOM *op);
void geom_apply(GEOM op, VECTOR *v);
float dvv(VECTOR v1, VECTOR v2);
float ddvv(VECTOR v1, VECTOR v2);
float avv(VECTOR v1, VECTOR v2);
float vdotv(VECTOR v1, VECTOR v2);
float avvvv(VECTOR v1, VECTOR v2, VECTOR v3, VECTOR v4);
VECTOR vector_normalize(VECTOR v);
VECTOR vector_vxv(VECTOR v1, VECTOR v2);
VECTOR vector_vplusv(VECTOR v1, VECTOR v2);
VECTOR vector_vminusv(VECTOR v1, VECTOR v2);
VECTOR vector_rescale(VECTOR v, float c);
VECTOR vector_sum3v(VECTOR v1, VECTOR v2, VECTOR v3);
VECTOR vector_neg(VECTOR v);
float dll(LINE line1, LINE line2);
float all(LINE line1, LINE line2);
float app(PLANE plane1, PLANE plane2);
LINE line_2v(VECTOR v1, VECTOR v2);
PLANE plane_3v(VECTOR v1, VECTOR v2, VECTOR v3);
float det3(float m[3][3]);
float det4(float m[4][4]);


/* File I/O functions */
STRINGS get_files(char *dir);
void free_strings(STRINGS *s);


/* Parameter functions */
int param_sav(char *key1, char *key2, char *key3, void *value, int s);
int db_close();
int db_open();
int iatom(char *conf_name, char *atom_name);
int param_get(char *key1, char *key2, char *key3, void *value);
int param_exist(char *key1, char *key2, char *key3);
int check_tpl(char *fname);
int load_param(char *fname);
int load_all_param(char *dirname);


#define LEN_KEY1 9
#define LEN_KEY2 6
#define LEN_KEY3 5
#define MAX_CONNECTED 8
#define MAX_CONNECTED2 MAX_CONNECTED *MAX_CONNECTED
#define MAX_CONNECTED3 MAX_CONNECTED2*MAX_CONNECTED

typedef struct {
    char ligand;
    char name[5];
    int  res_offset;
} CONNECTED_ATOM;

typedef struct {
    char orbital[10];
    int n;
    CONNECTED_ATOM atom[MAX_CONNECTED];
} CONNECT;

typedef struct {
    char   atom1[5];
    char   atom2[5];
    char   affected[MAXCHAR_LINE];
} ROTAMER;

typedef struct {
    int     n_swap;
    char   swap_atom1[8][5];
    char   swap_atom2[8][5];
} SWAP_RULE;    

typedef struct {
    char   atom1[5];
    char   atom2[5];
    char   atom3[5];
    int    n_term;
    float  V2[4];
    float  n_fold[4];
    float  gamma[4];
    int    opt_hyd;
} TORS;


/* data structure functions */
typedef struct {
    int    n;
    int    *res;
} MICROSTATE;

typedef struct {
   int ir;
   int ic;
   float value;
} INTERACTION;

typedef struct {
    int   sum_kc;
    int   sum_krkc;
    float E;
} UNIQ_STATE;

typedef struct ATOM_STRUCT {
    char   on;
    char   name[5];
    int    serial;
    char   resName[4];
    char   confName[6];
    char   chainID;
    int    resSeq;
    char   iCode;
    char   altLoc;
    char   confID[4];
    int    iConf;
    VECTOR xyz;
    float  crg;
    float  rad;
    float  sas;
    char   history[11]; /* 2 char for confName, Rotate, Swing, OH */
    struct ATOM_STRUCT *connect12[MAX_CONNECTED];
    int    connect12_res[MAX_CONNECTED];
    int    i_atom_prot;
    int    i_atom_conf;
    int    i_conf_res;
    int    i_res_prot;
    int    i_elem;
} ATOM;

typedef struct {
    char   on;
    char   uniqID[15];
    char   iCode;
    char   resName[4];
    int    resSeq;
    char   chainID;
    char   altLoc;
    char   confID[4];
    int    iConf;        /* Used to identify a conformer from mcce pdb,
                          * Also the index number to the head list
                          */
    char   confName[6];  /* conforname is resName+1st+2nd char in "history" */
    float  netcrg;       /* net charge */
    float  pKa;          /* pKa of the acid: HA <=> H+ + A- */
    float  Em;           /* Em of the reduced: RED <=> OX + e- */
    float  E_torsion;    /* torsion energy */
    float  E_vdw0;       /* VDW potential to its own (including backbone atoms) */
    float  E_vdw1;       /* VDW potential to backbone atoms (excluding its own) */
    float  E_epol;       /* elec to backbone */
    float  E_tors;
    float  E_rxn0;
    float  E_rxn;
    float  E_extra;      /* extra energy to make the calc. match training exep. */
    float  E_dsolv;
    float  E_ph;
    float  E_eh;
    float  E_mfe;
    float  E_self0;      /* self energy without mfe */
    double E_self;       /* total self energy, with mfe */
    float  occ;          /* occupance */
    int    e;            /* number of electron(s) gained */
    int    H;            /* number of proton(s) gained */
    char   history[10];  /* history of making this conformer */
    int    counter;      /* General counter */
    int    n_atom;
    ATOM   *atom;

    VECTOR r_min;
    VECTOR r_max;
    int    i_res_prot;
    int    i_conf_prot;  /* conformer index in the protein */
    int    i_conf_res;
    int    i_conf_subres;

    int    subresID;
    int    i_subres_res;

    float  tmp_pw_ele;
    float  tmp_pw_epol;
    float  tmp_pw_vdw;

    float  occ_old;
    char   toggle;
    float  *occ_table;
    int    counter_trial;
    int    counter_accept;
} CONF;

typedef struct {
    int    n_conf;
    CONF   **conf;

    int    k_subres;
} SUBRES;

typedef struct RES_STRUCT {
    int    resSeq;       /* used for identifying a residue */
    char   resName[4];   /* used for identifying a residue */
    char   chainID;      /* used for identifying a residue */
    char   iCode;        /* used for identifying a residue */
    int    cal_vdw;
    char   do_rot;       /* 0/1 flag to do rotamers */
    char   do_sw;
    int    opt_hyd;
    int    rotations;
    float  phi_swing;
    int    relax;
    int    n_hyd_pos;      /* number of position for each degree of freedom */
    float  sas;
    VECTOR r_min;
    VECTOR r_max;
    int    n_conf;
    CONF   *conf;
    int    n_conf_ori;    /* also used by jmao in step 3 to indicate the dielectric boundary conf */
    float  fixed_occ;     /* also used by jmao in step 3 to store pairwise interaction with boundary conf */

    /* used for step3 */
    int    i_bound;
    float  pw_bound;
    float  single_multi_bound;
    float  effective_epsilon;


    int **n_connect12;   /* (int) n_connect12[i_conf][i_atom] */
    int **n_connect13;
    int **n_connect14;
    ATOM ****connect12;  /* (ATOM *) connect12[i_conf][i_atom][i_connect] */
    ATOM ****connect13;
    ATOM ****connect14;
    float r12sq_max;
    float r13sq_max;
    float r14sq_max;

    CONF   *conf_w;
    CONF   *conf_new;
    CONF   *conf_old;

    int    i_res_prot;

    int    groupID;

    int    n_subres;
    SUBRES *subres;
    SUBRES *subres_w;
    SUBRES *subres_new;
    SUBRES *subres_old;
    int    i_subres_on;

    int    n_ngh;
    struct RES_STRUCT **ngh;

    int    nconf_limit;
    float  *sum_crg;
    int    n_flip_max;
    int    counter_trial;
} RES;

typedef struct {
    SUBRES  **subres;
} STATE;

typedef struct GROUP_STRUCT {
    int n_res;
    RES **res;

    int n_state;
    STATE *state;
    STATE *state_w;
    STATE *state_new;
    STATE *state_old;

    int n_ngh;
    struct GROUP_STRUCT **ngh;

    int n_flip_max;
} GROUP;

typedef struct {
    int n_res;
    RES *res;

    double E_base;
    double E_state;
    double E_min;
    double E_accum;
    double E_free_unfold;

    float  *H;
    float  *e;
    float  *sum_crg;

    int    n_group;
    GROUP  *group;
    int    nc;
    CONF   **conf;

    int  n_saved;
    UNIQ_STATE *saved_states;
} PROT;


typedef struct {
       float ele;
       float vdw;
       float crt;
       float ori;
       char  mark[4];
} PAIRWISE;


/*--- Data structure functions ---*/
ATOM   pdbline2atom(char *line);
PROT   load_pdb(FILE *fp);
int    write_pdb(FILE *stream, PROT prot);

int    ins_conf(RES *res, int ins, int n_atom);
int    del_conf(RES *res, int pos);
int    cpy_conf(CONF *tgt, CONF *src);

int    ins_res(PROT *prot, int ins);
int    del_res(PROT *prot, int pos);
int    cpy_res(RES *tgt, RES *src);

PROT   new_prot();
int    del_prot(PROT *prot);
int    cpy_prot(PROT *tgt, PROT *src);


/* Energy calculation */
float vdw(ATOM atom1, ATOM atom2);
VECTOR vdw_frc(VECTOR v1, VECTOR v2, float C6, float C12);
float vdw_conf(int ires, int iconf, int jres, int jconf, PROT prot);
float vdw_conf_hv(int ires, int iconf, int jres, int jconf, PROT prot);
float coulomb(ATOM atom1, ATOM atom2);
VECTOR coulomb_frc(VECTOR v1, VECTOR v2, float crg1, float crg2);
float coulomb_conf(int ires, int iconf, int jres, int jconf, PROT prot);
int setup_C6_C12(PROT prot);
void setup_vdw_fast(PROT prot);
void setup_connect_res(PROT prot, int i_res);
void free_connect_res(PROT prot, int i_res);
float vdw_conf_fast(int i_res, int i_conf, int j_res, int j_conf, PROT prot, int handle);
float coulomb_conf_fast(int i_res, int i_conf, int j_res, int j_conf, PROT prot);
int out_of_range(VECTOR i_min, VECTOR i_max, VECTOR j_min, VECTOR j_max, float range2);
float torsion_angle(VECTOR v0, VECTOR v1, VECTOR v2, VECTOR v3);
int torsion_atoms(CONF *conf_p, int i_atom, ATOM **atom0_p, ATOM **atom1_p, ATOM **atom2_p, ATOM **atom3_p, TORS *tors, int handle);
float torsion(float phi, float phi0, float n_fold, float barrier);
float torsion_conf(CONF *conf_p);
VECTOR torsion_torq(float phi, float V2, float n_fold, float gamma, VECTOR k);
float Ecoulomb_conf2conf(PROT prot, int ir, int ic, int kr, int kc, float epsilon);
float Evdw_conf2conf(PROT prot, int ir, int ic, int kr, int kc);
float CoulombBySAS(ATOM atom1, ATOM atom2);
float RxnBySAS(CONF conf);
int sas_native(PROT prot);

/* other functions */
int strip(char *target, char *str);
int get_env();
int assign_rad(PROT prot);
int assign_crg(PROT prot);
int get_connect12(PROT prot);
int get_connect12_conf(int i_res, int i_conf, PROT prot);
int surfw(PROT prot, float probe_rad);
int surfw_res(PROT prot, int ir, float probe_rad);
int surfw_l2(PROT prot, float probe_rad);
void shuffle_n(int *array, int n);
int cmp_conf(CONF conf1, CONF conf2, float IDEN_THR);
int cmp_conf_hv(CONF conf1, CONF conf2, float IDEN_THR);
void id_conf(PROT prot);
int sort_conf(PROT prot);
int sort_res(PROT prot);
void get_vdw0(PROT prot);
void get_vdw1(PROT prot);
void load_headlst(PROT prot);
void write_headlst(PROT prot);
double ran2(long *idum);
int cmp_Eself(const void *a, const void *b);
float hbond_extra(CONF a, CONF b);
int def(FILE *source, FILE *dest, int level);
int inf(FILE *source, FILE *dest);
int del_dir(char *dir_name);
int make_matrices(PROT prot, char *dir);
int load_energies(PAIRWISE **pairwise_raw, int n);
int add_matrices(PAIRWISE **pw1, PAIRWISE **pw2);

/* Modeules */
int init();
int premcce();
int rotamers();
int energies();
int energies2();
int monte();
int monte2();

