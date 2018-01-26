function get_gene_info()
global sample_base_dir;
global sample_num;
global gene_body_len;
global gene_id;
global gene_name;
global cancer_name;
global cancer_stage;
global figdir;
global tss_len;
global out_methy_file;
global out_ins_file;
global out_del_file;
global out_snp_file;
global out_ins_num_file;
global out_del_num_file;
global out_snp_num_file;
global out_ave_methy_file;

%gene_id = 831; %APC
gene_id = 5572; %FAT4
gene_name = 'FAT4';
tss_len = 2000;
cancer_name='COAD';
cancer_stage='i';
base_path = '../../figures/';

%file_name_definition
out_methy_file = 'out_methy.txt';
out_ins_file = 'out_INS.txt';
out_del_file = 'out_DEL.txt';
out_snp_file = 'out_SNP.txt';
out_ins_num_file = 'out_INS_num.txt';
out_del_num_file = 'out_DEL_num.txt';
out_snp_num_file = 'out_SNP_num.txt';
out_ave_methy_file = 'out_ave_methy.txt';

%count gene_body_len
gene_label_path = '../../global_files/gene_label.dat';
gene_label_info = load(gene_label_path);
gene_body_absolute_start = gene_label_info(gene_id, 9);
gene_body_absolute_end = gene_label_info(gene_id, 10);
gene_body_len = gene_body_absolute_end - gene_body_absolute_start;

py.importlib.import_module('get_gene_data');

%count_sample_num
sample_base_dir = 'G:/intermediate_file/common_patients_data/';
out_sample_num_file = 'out_sample_num.txt';
py.get_gene_data.get_sample_num(sample_base_dir, cancer_name, cancer_stage);
sample_num = load(out_sample_num_file);

figdir = strcat(base_path, 'mutation_classification/', cancer_name,'/');
if ~exist(figdir)
    mkdir(figdir);
end

end