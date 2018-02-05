function get_gene_info(cancer_name, gene_name, gene_id)
global sample_base_dir;
global sample_num;
global out_sample_num_file;
global normal_sample_num;
global gene_body_len;
global cancer_stage;
global figdir;
global tss_len;
global out_methy_file;
global out_normal_methy_file;
global out_ins_file;
global out_del_file;
global out_snp_file;
global out_ins_num_file;
global out_del_num_file;
global out_snp_num_file;
global out_ave_methy_file;
global normal_mean_methy_file;
global normal_stage;
tss_len = 2000;
cancer_stage='i';
normal_stage = 'normal';
base_path = '../../figures/';

%file_name_definition
out_methy_file = 'out_methy.txt';
out_normal_methy_file = 'out_normal_methy.txt';
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

Path = py.sys.path;
if count(Path,'.') == 0
    insert(Path,int32(0),'.');
end

py.importlib.import_module('get_gene_data');

%count_sample_num
% sample_base_dir = 'G:/intermediate_file/common_patients_data/';
sample_base_dir = '/Volumes/Elements/intermediate_file/common_patients_data/';
normal_mean_methy_file = strcat(sample_base_dir, cancer_name, '/', normal_stage, '/normal_mean_methy/', gene_name, '.tsv');

out_sample_num_file = 'out_sample_num.txt';
py.get_gene_data.get_sample_num(sample_base_dir, cancer_name, cancer_stage, out_sample_num_file);
sample_num = load(out_sample_num_file);

out_normal_sample_num_file = 'out_normal_sample_num.txt';
py.get_gene_data.get_sample_num(sample_base_dir, cancer_name, normal_stage, out_normal_sample_num_file);
normal_sample_num = load(out_normal_sample_num_file);

figdir = strcat(base_path, 'mutation_classification/', cancer_name, '/', gene_name, '/');
if ~exist(figdir)
    mkdir(figdir);
end

end