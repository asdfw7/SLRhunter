import sys,os
import copy
import argparse
import shutil
import json
import math
import multiprocessing
from collections import OrderedDict, Counter
from xopen import xopen as open
from Bio import SeqIO
#from . import Seqs
from .KMC import run_kmc_dict, KmerMatrix
#from .Cluster import Cluster
#from . import Stats
#from . import Circos
#from . import Blocks
from .small_tools import mkdirs, rmdirs, mk_ckp, check_ckp, test_s
from .RunCmdsMP import logger, available_memory, limit_memory
from .__version__ import version


bindir = os.path.dirname(os.path.realpath(__file__))
NCPU = multiprocessing.cpu_count()
MEM = available_memory()

def makeArgparse():
	parser = argparse.ArgumentParser( 
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description='Detect SDR based on kmers.',
		)
	# input
	group_in = parser.add_argument_group('Input', 'Input genome and config files')
	group_in.add_argument('-s',  dest='sd', required=True, metavar='FILE',
					help="sample design file [required]")
	group_in.add_argument('-g',  dest='genome', default=None, metavar='FILE',
					help="genome fasta file")
	group_in.add_argument('-subset',  metavar='FILE',
					help="subset of sample design file")

	# output
	group_out = parser.add_argument_group('Output')
	group_out.add_argument('-pre', '-prefix', default=None, dest='prefix', metavar='STR',
					help="Prefix for output [default=%(default)s]")
	group_out.add_argument('-o', '-outdir', default='sdr-results', dest='outdir', metavar='DIR',
					help="Output directory [default=%(default)s]")
	group_out.add_argument('-tmpdir', default='tmp', type=str, metavar='DIR',
					help="Temporary directory [default=%(default)s]")
	# kmer
	group_kmer = parser.add_argument_group('Kmer', 'Options to count and filter kmers')
	group_kmer.add_argument('-k', type=int, default=31, metavar='INT',
					 help="Length of kmer [default=%(default)s]")
	group_kmer.add_argument('-f', '-min_fold', type=float, default=1.5, metavar='FLOAT', dest='min_fold',
					help="Minimum fold [default=%(default)s]")
	group_kmer.add_argument('-baseline', type=int, default=1, 
					 help="Use sub-maximum (1) or minimum (-1) as the baseline of fold \
[default=%(default)s]")

	group_kmer.add_argument('-lower_count', type=int, default=3, metavar='INT',
					 help="Don't output k-mer with count < lower-count [default=%(default)s]")
	group_kmer.add_argument('-low_mem', action="store_true", default=None,
					help="Low MEMory but slower [default: True if genome size > 3G, else False]")
	group_kmer.add_argument('-by_count', action="store_true", default=False,
					help="Calculate fold by count instead of by proportion [default=%(default)s]")
	group_kmer.add_argument('-mapq', type=int, default=0, metavar='INT',
                     help="Minimum mapping quality of kmers, 1 to remove multiple hits [default=%(default)s]")

	# filter
 
	# cluster
	group_clst = parser.add_argument_group('Cluster', 'Options for clustering to phase')
	group_clst.add_argument('-max_pval', type=float, default=0.05, metavar='FLOAT',
					help="Maximum P value for all hypothesis tests [default=%(default)s]")
	group_clst.add_argument("-test_method", default='kruskal', 
					choices=['ttest_ind', 'kruskal', 'wilcoxon', 'mannwhitneyu'],
					help="The test method to identify differiential kmers[default=%(default)s]")
					
	group_clst.add_argument("-figfmt", default='pdf', type=str, 
					choices=['pdf', 'png', ], # 'svg','tiff', 'jpeg', 'bmp'],
					help="Format of figures [default=%(default)s]")
	group_clst.add_argument('-heatmap_colors', nargs='+', default=('blue', 'yellow', 'red'), metavar='COLOR',
					help="Color panel (2 or 3 colors) for heatmap plot [default: %(default)s]")
	group_clst.add_argument('-heatmap_options', metavar='STR',
					default="Rowv=T,Colv=T,scale='col',dendrogram='row',labCol=F,trace='none',\
key=T,key.title=NA,density.info='density',main=NA,xlab='Differential kmers',margins=c(2.5,12)",
					help='Options for heatmap plot (see more in R shell with `?heatmap.2` \
of `gplots` package) [default="%(default)s"]')


	# circos
	group_circ = parser.add_argument_group('Circos', 'Options for circos plot')
	group_circ.add_argument('-disable_circos', action="store_true", default=True,
					help="Disable this step [default=%(default)s]")
	group_circ.add_argument('-window_size', type=int, default=1000000, metavar='INT',
					help="Window size (bp) for circos plot [default=%(default)s]")
	group_circ.add_argument('-disable_blocks', action="store_true", default=False,
					help="Disable to plot homologous blocks [default=%(default)s]")
	group_circ.add_argument("-aligner", metavar='PROG', 
					default='minimap2', 
					choices=['minimap2', 'unimap'],
					help="Programs to identify homologous blocks [default=%(default)s]")
	group_circ.add_argument("-aligner_options", metavar='STR',
					default='-x asm20 -n 10',
					help='Options for `-aligner` to align chromosome sequences [default="%(default)s"]')
	group_circ.add_argument('-min_block', type=int, default=100000, metavar='INT',
					help="Minimum block size (bp) to show [default=%(default)s]")
	group_circ.add_argument('-alt_cfgs', nargs='+', metavar='CFGFILE', default=None,
					help="An alternative config file for identifying homologous blocks \
[default=%(default)s]")
	group_circ.add_argument("-chr_ordered", default=None, type=str, metavar='FILE',
					help="Provide a chromosome order to plot circos \
[default=%(default)s]")

	# others
	group_other = parser.add_argument_group('Other options')
	group_other.add_argument('-p', '-ncpu', type=int, default=NCPU, metavar='INT', dest='ncpu',
					 help="Maximum number of processors to use [default=%(default)s]")
	group_other.add_argument('-max_memory', type=str, default=MEM, metavar='MEM', 
					 help="Maximum memory to use where limiting can be enabled. [default=%(default)s]")
	group_other.add_argument('-cleanup', action="store_true", default=False,
					help="Remove the temporary directory [default=%(default)s]")	
	group_other.add_argument('-overwrite', action="store_true", default=False,
					help="Overwrite even if check point files existed [default=%(default)s]")
	group_other.add_argument('-v', '-version', action='version', version=version)
	
	args = parser.parse_args()
	if args.prefix is not None:
		args.prefix = args.prefix.replace('/', '_')
		args.outdir = args.prefix + args.outdir
		args.tmpdir = args.prefix + args.tmpdir
	return args

		
class Pipeline:
	def __init__(self, sd, **kargs):
		self.sdfile = sd
		self.sd = Design(sd)
		self.grouped = self.sd.groupby_group()	#{'g1': ['s1', 's2'], 'g2': ['s3', 's4']}
		self.groups = self.grouped.keys()	# ['g1', 'g2']
		self.grouped_samples = self.grouped.values()	# [['s1', 's2'], ['s3', 's4']]
		self.samples =  self.sd.samples			# ['s1', 's2', 's3', 's4']
		self.group_dict = self.sd.group_dict	# {'s1': 'g1', 's2': 'g1', 's3': 'g2', 's4': 'g2'}
		self.nsamples = len(self.samples)
		self.__dict__.update(**kargs)
		if self.subset:
			self.subset = Design(self.subset).samples
			
		self.kargs = kargs

		
	def run(self):
		self.pool_method = 'imap_unordered'
		# mkdir
		self.outdir = os.path.realpath(self.outdir)
		self.tmpdir = os.path.realpath(self.tmpdir)
		mkdirs(self.outdir)
		mkdirs(self.tmpdir)
		self.outdir += '/'
		self.tmpdir += '/'
		self._outdir = self.outdir
		if self.prefix is not None:
			self.outdir = self.outdir + self.prefix
			self.tmpdir = self.tmpdir + self.prefix

		logger.info('###Step: Kmer Count')
		# jellyfish
		outdir = '{}/kmer'.format(self.tmpdir)
		mkdirs(outdir)
		logger.info('Counting kmer by KMC')
		merged_dumpfile, total_kmers = run_kmc_dict(self.sd.datafile_dict, 
						k=self.k, lower_count=self.lower_count,
						outdir=outdir, threads=self.ncpu, overwrite=self.overwrite)	# ckp for each sample
		logger.info('Total kmers: {}'.format(total_kmers))

		# matrix
		logger.info('Loading kmer matrix from jellyfish')	# multiprocessing by kmer
		chunksize = None if self.pool_method == 'map' else 1000
#		dumps = JellyfishDumpSamples(merged_dumpfile, samples=self.samples, groups=self.grouped,
#							ncpu=self.ncpu, )
							#method=self.pool_method, chunksize=chunksize)
		self.basename = 'k{}'.format(self.k, )
		self.para_prefix = '{}{}'.format(self.outdir, self.basename)
#		sig_matfile = self.para_prefix + '.sig.kmer.mat'
		matfile = self.para_prefix + '.kmer.mat'
#		sig_kmers = self.para_prefix + '.sig.kmer-group.tsv'

		ckp_file = self.mk_ckpfile(matfile)
		ckp = check_ckp(ckp_file)
		if self.overwrite or not ckp or not test_s(matfile):
			mat = KmerMatrix(merged_dumpfile, ncpu=self.ncpu,method=self.pool_method, chunksize=chunksize)
			with open(matfile, 'w') as fout:
				mat.filter(fout, 
					subset=self.subset, d_groups=self.grouped, 
					min_mean_freq=1, max_missing_rate=0.3, by_group=True, 
					min_maf=0.1)
			mk_ckp(ckp_file, self.subset)
		#else:
		#	d_kmers, = ckp
		mat = KmerMatrix(matfile, ncpu=self.ncpu,method=self.pool_method, chunksize=chunksize)
		logger.info('Generating data for GEMMA')
		mat.to_gemma(matfile, tmpdir=self.tmpdir, genome=self.genome, mapq=self.mapq)
		
		return
		# kmeans cluster
		# logger.info('###Step: Cluster')
		# cluster = Cluster(matfile, sg_assigned=self.sd.group_dict,
				# re_assign=False, bootstrap=False)
		# self.d_sg = d_sg = cluster.d_sg	# chrom -> SG
		# logger.info('Subgenome assignments: {}'.format(d_sg))

		# specific kmers and location
		# sg_kmers = self.para_prefix + '.sig.kmer-subgenome.tsv'
		# logger.info('Outputing significant differiential `kmer` - `subgenome` maps to `{}`'.format(sg_kmers))
		# with open(sg_kmers, 'w') as fout:	# multiprocessing by kmer
			# # kmer -> SG
			# d_kmers = cluster.output_kmers(fout, max_pval=self.max_pval, ncpu=self.ncpu,
							# test_method=self.test_method)
		# logger.info('{} significant specific kmers'.format(len(d_kmers)//2))
		# for sg, count in sorted(Counter(d_kmers.values()).items()):
			# logger.info('\t{} {}-specific kmers'.format(count//2, sg))
		
		
		# heatmap	# static method
		outfig = dumps.heatmap(sig_matfile, mapfile=self.sdfile, kmermapfile=sig_kmers,
					figfmt=self.figfmt, color=self.heatmap_colors, 
					heatmap_options=self.heatmap_options)
		
		# PCA
		# outfig = self.para_prefix + '.kmer_pca.' + self.figfmt
		# logger.info('Outputing PCA plot to `{}`'.format(outfig))
		# cluster.pca(outfig, n_components=self.nsg)

		# return

		# if self.just_core:
			# self.step_final()
			# logger.info('Pipeline completed early')
			# return
		
		sg_map = self.para_prefix + '.group.bin.count'
		ckp_file = self.mk_ckpfile(sg_map)
		if self.overwrite or self.re_filter or not check_ckp(ckp_file) or not test_s(sg_map):	
		# SG id should be stable
			logger.info('Outputing `coordinate` - `group` maps to `{}`'.format(sg_map))
			chunksize = None if self.pool_method == 'map' else 10
			with open(sg_map, 'w') as fout:	# multiprocessing by chrom chunk
				Seqs.map_kmer3(chromfiles, d_kmers, fout=fout, k=self.k, 
					bin_size=10000, sg_names=self.sg_names,
					ncpu=self.ncpu, method=self.pool_method, chunksize=chunksize)
			mk_ckp(ckp_file)
		# enrich by BIN
		logger.info('Enriching group by chromosome window (size: {})'.format(self.window_size))
	#	bins, counts = Circos.counts2matrix(sg_map, keys=self.sg_names, keycol=3, window_size=self.window_size)
		bins, counts = Circos.stack_matrix(sg_map, window_size=self.window_size)
	#	logger.info('Matrix loaded')
		bin_enrich = self.para_prefix + '.bin.enrich'
		bin_exchange = self.para_prefix + '.bin.group'
		with open(bin_enrich, 'w') as fout:	# multiprocessing by chrom bin
			with open(bin_exchange, 'w') as fout2:	# group exchanges
				sg_lines = Stats.enrich_bin(fout, fout2, self.d_sg, counts, colnames=self.sg_names, rownames=bins,
					max_pval=self.max_pval, ncpu=self.ncpu)
		logger.info('Output: {}'.format(bin_enrich))

		# # custom
		# if self.custom_features is not None:
	# #		for i, feature in enumerate(self.custom_features):
			# feat_map = self.para_prefix + '.custom.bin.count'
			# logger.info('Mapping subgenome-specific kmers to custom features: {}'.format(self.custom_features))
			# pool_method = self.pool_method
			# chunksize = None if pool_method == 'map' else 10000
			# with open(feat_map, 'w') as fout:	# multiprocessing by LTR
				# Seqs.map_kmer3(self.custom_features, d_kmers, fout=fout, k=self.k, ncpu=self.ncpu, 
								# bin_size=10000000, sg_names=self.sg_names,
								# chunk=False, log=False, method=pool_method, chunksize=chunksize)
			# # enrich SG by feature
			# logger.info('Enriching subgenome-specific features')
			# bins, counts = Circos.stack_matrix(feat_map, window_size=100000000)
			# feat_enrich = self.para_prefix + '.custom.enrich'
			# with open(feat_enrich, 'w') as fout:
				# d_enriched, *_ = Stats.enrich_ltr(fout, self.d_sg, counts, colnames=self.sg_names, rownames=bins, 
						# max_pval=self.max_pval, ncpu=self.ncpu)
			# logger.info('Output: {}'.format(feat_enrich))
			
			# logger.info('{} significant subgenome-specific features'.format(len(d_enriched)))
			# for sg, count in sorted(Counter(d_enriched.values()).items()):
				# suffix = {'shared':''}.get(sg, '-specific')
				# logger.info('\t{} {}{} features'.format(count, sg, suffix))
		
		# # LTR
		# ltr_bedlines, enrich_ltr_bedlines = self.step_ltr(d_kmers) if not self.disable_ltr else ([],[])

		# circos
		if self.chr_ordered:
			self.chromfiles = [self.d_chromfiles[chr] for chr in self.chr_ordered]
		if not self.disable_circos:
			self.step_circos(
				bedfile=sg_map, # chrom - coord - SG			circles 3 - n+2
				sg_lines=sg_lines, # SG ratio and enrichment	circles 1 - 2
				d_sg = self.group_dict, # sample -> SG, for colors
				prefix=self.para_prefix + '.circos',
				figfmt=self.figfmt,
				window_size=self.window_size)

		self.step_final()
		logger.info('Pipeline completed')

	def mk_ckpfile(self, file):
		return '{}{}.ok'.format(self.tmpdir, os.path.basename(file))


	def step_circos(self, *args, **kargs):
		logger.info('###Step: Circos')
		# blocks
		if not self.disable_blocks:
			pafs, paf_offsets = self.step_blocks()
			kargs['pafs'] = pafs
			kargs['paf_offsets'] = paf_offsets
			kargs['min_block'] = self.min_block

		# circos
		circos_dir = bindir+'/circos'
		wkdir = self.para_prefix + '.circos' #self._outdir+'/circos'
		rmdirs(wkdir)
		try:
			logger.info('Copy `{}` to `{}`'.format(circos_dir, self._outdir))
			shutil.copytree(circos_dir, wkdir)
		except FileExistsError:
			pass
		Circos.circos_plot(self.chromfiles, wkdir, *args, **kargs)
	
	def step_blocks(self):
		outdir = '{}Blocks/'.format(self.tmpdir)
		multiple = {'minimap2': 20, 'unimap': 40}	# relative with repeat amount
		max_size = max(self.d_size.values())
		mem_per_cmd = max_size*math.log(max_size, 10) * multiple[self.aligner]
		ncpu = min(self.ncpu, limit_memory(mem_per_cmd, self.max_memory))
		logger.info('Using {} processes to align chromosome sequences'.format(ncpu))
		thread = int(self.ncpu // ncpu)
		
		mkdirs(outdir)
		pafs, paf_offsets = Blocks.run_align(self.alt_sgs, self.d_chromfiles, outdir, aligner=self.aligner,
						ncpu=ncpu, thread=thread, d_size=self.d_size, overlap=self.min_block*5,
						opts=self.aligner_options, overwrite=self.overwrite)
		#print(paf_offsets)
		return pafs, paf_offsets
	
	def step_final(self):
		# cleanup
		if self.cleanup:
			logger.info('Cleaning {}'.format(self.tmpdir))
			rmdirs(self.tmpdir)
	def sort_labels(self, order, labels, chromfiles):
		d = dict(zip(labels, chromfiles))
		labels, chromfiles = [], []
		for lab in order:
			if not lab in d:
				continue
			chromfile = d[lab]
			labels += [lab]
			chromfiles += [chromfile]
		return labels, chromfiles
		
class Design:
	def __init__(self, design=None, samples=None):
		self.design = design
		if samples is None and design is not None:
			self._samples = list(self._parse())
		else:
			self._samples = samples
	@property
	def samples(self):
		return [sample.sample for sample in self._samples]
	@property
	def groups(self):
		return [sample.group for sample in self._samples]
	@property
	def datafiles(self):
		return [sample.datafiles for sample in self._samples]
	@property
	def sample_dict(self):
		return OrderedDict((sample.sample, sample) for sample in self._samples)
	@property
	def group_dict(self):
		return OrderedDict((sample.sample, sample.group) for sample in self._samples)
	@property
	def datafile_dict(self):
		return OrderedDict((sample.sample, sample.datafiles) for sample in self._samples)
	def groupby_group(self):
		d_groups = OrderedDict()
		for sample in self._samples:
			group = sample.group
			sample = sample.sample
			try: d_groups[group].append(sample)
			except KeyError: d_groups[group] = [sample]
		return d_groups
	def __iter__(self):
		return iter(self._samples)
	def __len__(self):
		return len(self._samples)
	def append(self, sample):
		self._samples.append(sample)
	def _parse(self):
		for line in open(self.design):
			yield DesignLine(line)
class DesignLine:
	def __init__(self, line):
		temp = line.strip().split()
		self.sample, self.group = temp[:2]
		self.datafiles = temp[2:]
	def __str__(self):
		return self.sample

def parse_idmap(mapfile=None):
	'''idmap: old_id new_id'''
	if not mapfile:
		return None
	d_map = OrderedDict()
	for line in open(mapfile):
		line = line.strip().split('#')[0]
		if not line:
			continue
		temp = line.split()
		old_id = temp[0]
		try: new_id = temp[1]
		except IndexError: new_id = old_id.split('|')[-1]
		d_map[old_id] = new_id
	return d_map
def check_duplicates(lst):
	count = Counter(lst)
	duplicates = {v: c for v,c in count.items() if c>1}
	if duplicates:
		raise ValueError('Duplicates detected: {}'.format(duplicates))

class SGConfig:
	def __init__(self, sgcfg, **kargs):
		self.sgcfg = sgcfg
		self.kargs = kargs
		self.sgs = list(self)
	def __iter__(self):
		return self._parse()
	def _parse(self):
		self.nsgs = []
		self.nsg = 0
		self.chrs = []
		for line in open(self.sgcfg):
			temp = line.split('#')[0].strip().split()
			if not temp:
				continue
			chrs = [list(map(lambda x: add_prefix(x, **self.kargs), v.strip(',').split(','))) 
						for v in temp]
			self.nsgs += [len(chrs)]
			if self.nsg == 0:
				self.nsg = len(chrs)
			if len(chrs) != self.nsg:
				logger.warn('Number of column is different in line {}: \
{} in this line but {} in previous line'.format(temp, len(chrs), self.nsg))
			for xchr in chrs:
				for xxchr in xchr:
					self.chrs += [xxchr]
			yield chrs
		self.nsg = max(self.nsgs)
		for chr, count in Counter(self.chrs).items():
			if count > 1:
				logger.warn('Chromsome id {} repeat {} times'.format(chr, count))

def add_prefix(val, prefix=None, sep='|'):
	if prefix:
		vals = ['{}{}'.format(prefix, v) for v in val.split(sep) if v]
		return ''.join(vals)
	else:
		return val

def main():
	args = makeArgparse()
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	pipeline = Pipeline(**args.__dict__)
	pipeline.run()


if __name__ == '__main__':
	main()

