#coding: utf-8

import requests

from bs4 import BeautifulSoup  #解析网页，获取标签下面的内容
import json
import lxml
import re
import time
import random

base_url = 'https://www.ncbi.nlm.nih.gov'
gene_search_url = 'https://www.ncbi.nlm.nih.gov/gene/?term='

headers = {
	'Accept':'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,image/apng,*/*;q=0.8',
	'Accept-Encoding':'gzip, deflate, br',
	'Accept-Language':'zh-CN,zh;q=0.9',
	'Connection':'keep-alive',
	'Cookie':'_gat=1; starnext=MYGwlsDWB2CmAeAXAXAJgLy2ogTrAXgGQDM6A5lrIQCzoAWA9gLZUCMtADgIYUD6AbmFgB3QqwDs6AGZcQAZzYAOabIWFUABnQBhAKKLdqRQFZiAQQBsugELaAnMdYbnL1wDFNGi3YB0TXuLqrJjYePgApOYUcJFmjCyxAHIA8om66hjCWT7QwABGYDkgTDlgdD5kDPzqtFi4BLHRVKgWIfURUZSx8c2SdWGNXeY9SanpRm0D4tpN4dMj0ylp6nborHaaJFpOihokweubxBhO1HvEpDLyVMS0uACuN8ZrxqioJK2orIrUJJKov2IyicxnEGnEFgsxhIqxBYNQxmoxBo2w0oJ2XxowXEiMUNAw0C4iDA/Co1FooAgMAQiBoz3e1FaVzU1EkrFYxmh1GUiJoqzy90QiAY0F4Ci4OGAdHCqG00AYwhwXA4hEc6HFkroqowAGVYBKparSBziBZVbRTHZCFxcrA5MKcABJAAm6F0oQIbgYOCY1tt9u92hAXDkckSXBY6HliuVABoNVLeFJvUx4/rNbG6PrnbAcABuDh4fgJuhoSYEEjkSg0ejMNicHiwARCUQSFTXMTKZnNLR6AxGUyWGz2RyuMcaDzObx+AJBcsdMyzYZ10bLN7oLLCHL5QrQYqlcqVaoA+eDODqVr9BqdGLLljqPoehdLuIr8xLcbKK8RaazeYrxYxhWNYNnOVFdn2ECjhODQzhIS5VBuO4cEeEhng5N4PnQL4fj+bDAWBNEwQhXDiFhIiNARCxkTONYiJ2ahWCxdAcW5fF0GVMAaFoEs4CQOl0GRUwXlNc0WPEK0LAwJpCGo9AzAABW0PN7j3BguGdPiUAwb9Kxk2gejEDBVJAdTnTEBs+EEEQxEkbtO3bNRNB0fRDBMcwrFsBwnHHFxJy8Xx/ECL5TxvWBujfMwPwyDdslyAoihKaAygqKoalCxdqxaDKX0M1BH3aM9wrvYqoqAiZvzmGYun/BJALXWFQK2OiIOIA4muOOi4IuRykPQB4nheTDTWw75AX+Ai6NBcFvDxMipvhUFAlouEdk5ZjWLxagCSJEkyVoEyzK0gToWEk0zURcSrU5dBFAu1odj2UFsKe5QcWMa7Vj2CwtG+9riHEIA==',
	'Host':'www.ncbi.nlm.nih.gov',
	'Referer':'https://www.ncbi.nlm.nih.gov/gene',
	'Upgrade-Insecure-Requests':'1',
	'User-Agent':'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/62.0.3202.94 Safari/537.36'
}

gene_name_list = ['BRAF', 'FLG', 'KRAS', 'MTOR', 'MUC4', 'PIK3CA', 'AHNAK', 'APC', 'ATM', 'CREBBP', 'CSMD1', 'DLC1', \
	'DMD', 'FAT1', 'FAT4', 'LRP1B', 'MYH9', 'PCDH17', 'SMAD4', 'TP53', 'ZFHX3', 'GLI1', 'ALK', 'CTNND2', 'DLX5', 'EGFR', \
	'FEZF1', 'FGF4', 'FGF5', 'LMO1', 'OTX2', 'PRDM14', 'TAL1', 'AKR1B1', 'BTG4', 'CHST10', 'CNTNAP2', 'DCC', \
	'NRCAM', 'PTPRD', 'PTPRT', 'SFRP2', 'THBD', 'TMEFF2', 'UNC5A']

#获取主页源码
def get_gene_info_link(origin_url, gene_name):
	url = origin_url + gene_name
	html = requests.get(url, headers=headers).text  #发送请求获取主页文本
	#获取到的html为str类型，已经是unicode类型了
	html = html.encode('utf-8')
	#print html
	soup = BeautifulSoup(html,'lxml')  #解析网页，第二个参数为解析方式，'html.parse'为python自带，速度较慢，lxml直接解析网页，速度更快
	links = [link.find('td', class_='gene-name-id').find('a').get('href') for link in soup.find_all('tr', class_='rprt') if (link.find('em').text == u'Homo sapiens')]
	gene_sub_url = links[0]  #links里面会存放人类的gene对应的href
	print gene_sub_url
	return base_url + gene_sub_url

def get_gene_summary(gene_info_url):
	html = requests.get(gene_info_url, headers=headers).text  #发送请求获取主页文本
	
	#获取到的html为str类型，已经是unicode类型了
	html = html.encode('utf-8')
	soup = BeautifulSoup(html,'lxml')
	summary = [link.find_next_sibling().string for link in soup.find_all('dt') if (link.string == u'Summary')]
	if summary == None or len(summary) == 0:
		summary = ['none']
	return summary[0]

def get_gene_info(gene_name, gene_summary_list=None):
	gene_info_url = get_gene_info_link(gene_search_url, gene_name)
	summary = get_gene_summary(gene_info_url)
	print summary
	if(gene_summary_list != None):
		gene_summary_list.append(summary)
	result_pro = re.findall(r'proliferation', summary)
	result_dif = re.findall(r'differentiation', summary)
	if len(result_pro) > 0 and len(result_dif) > 0:
		return 'both'
	elif len(result_pro) > 0:
		return 'proliferation'
	elif len(result_dif) > 0:
		return 'differentiation'
	else:
		return 'neither'

def write_gene_list_summary(file_name, gene_summary_list):
	with open(file_name, 'w') as f:
		for idx in range(len(gene_summary_list)):
			gene_name = gene_name_list[idx]
			gene_summary = gene_summary_list[idx]
			line_str = gene_name + '\t' + gene_summary + '\n'
			f.write(line_str)

def write_gene_category_file(file_name, gene_list):
	with open(file_name, 'w') as f:
		for gene_name in gene_list:
			line_str = gene_name + '\n'
			f.write(line_str)

def write_gene_category_list(gene_info_list):
	gene_info_set = set(gene_info_list)
	gene_info_dict = {}
	for idx in range(len(gene_info_list)):
		gene_name = gene_name_list[idx]
		gene_info = gene_info_list[idx]
		gene_info_dict[gene_info] = gene_info_dict.get(gene_info, [])
		gene_info_dict[gene_info].append(gene_name)
	print gene_info_dict
	for gene_category in gene_info_dict:
		file_name = gene_category + '.txt'
		sub_name_list = gene_info_dict[gene_category]
		write_gene_category_file(file_name, sub_name_list)


def get_gene_info_list(summary=False):
	gene_info_list = []
	gene_summary_list = None
	if(summary):
		gene_summary_list = []
	for gene_name in gene_name_list:
		gene_info = get_gene_info(gene_name, gene_summary_list)
		time.sleep(random.randint(1, 5))
		gene_info_list.append(gene_info)
	if(summary):
		return [gene_info_list, gene_summary_list]
	return gene_info_list

