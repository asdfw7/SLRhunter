#!/bin/env python3
# py3
import sys
import numpy as np
from collections import Counter
class Gemma:
	def __init__(self, output):
		self.output = output
	def __iter__(self):
		return self._parse()
	def _parse(self):
		i = 0
		for line in open(self.output):
			i += 1
			if i == 1:
				title = line
				continue
			yield GemmaLine(line, title)
class GemmaLine:
	def __init__(self, line, title):
		self.title = title.strip()
		self.line = line.strip()
		self._parse()
	def _parse(self):
		title = self.title.split()
		line = self.line.split()
		for k, v in zip(title, line):
			setattr(self, k, tr_numeric(v))

class BamCov:
	def __init__(self, incov, chroms=[], maf=0.05,  **parse_arg):
		self.incov = incov
		self.chroms = set(chroms)
		self.maf = maf
		self.parse_arg = parse_arg
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.incov):
			if not line.strip():
				continue
			cov = BamCovLine(line, **self.parse_arg)
			if cov.maf < self.maf:
				continue
			yield cov
	def to_matrix(self):
		coord = []
		data = []
		for rc in self:
			if self.chroms and rc.chrom not in self.chroms:
				continue
			coord += [(rc.chrom, rc.start, rc.end)]
			data += [rc.covs]
		return coord, np.array(data), 
class BamCovLine:
	def __init__(self, line, noend=False):
		chrom, start, end, *covs = line.strip().split()
		self.chrom = chrom
		self.start = tr_numeric(start)
		self.end = tr_numeric(end)
		self.covs = list(map(tr_numeric, covs))
		if noend:
			self.covs = [self.end] +self.covs
			self.end = None
	def __len__(self):
		return self.end - self.start
	def __str__(self):
		return '{chrom}:{start}-{end}'.format(**self.__dict__)
	@property
	def maf(self):
		return 1 - 1.0*len([cov for cov in self.covs if cov==0])/len(self.covs)
class BimBam:
	def __init__(self, bbfile):
		self.bbfile = bbfile
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.bbfile):
			yield BBLine(line)
	def to_matrix(self):
		data = []
		for rc in self:
			data += [rc.values]
		return np.array(data)
	def get_column(self, n=1):
		data = self.to_matrix()
		return data[:, n-1]

class BBLine:
	def __init__(self, line, sep='\t'):
		values = line.strip().split(sep)
		self.values = list(map(tr_numeric, values))	

def tr_numeric(val):
	try: return int(val)
	except:
		try: return float(val)
		except: return val
def reorder_cov(incov, old_sd, new_sd, outcov):
	old_samples = [line.strip().split()[0] for line in open(old_sd)]
	new_samples = [line.strip().split()[0] for line in open(new_sd)]
	new_idx = [old_samples.index(s) for s in new_samples]
	print(new_idx, file=sys.stderr)
	coord, data = BamCov(incov).to_matrix()
	#print(data, file=sys.stderr)
	data = data[:, new_idx]
	#print(data, file=sys.stderr)
	for (chrom, start, end), covs in zip(coord, data):
		line = [chrom, start, end] + list(covs)
		line = map(str, line)
		line = '\t'.join(line)
		print(line, file=outcov)

def cov2gemma(incov, outgeno, fold=2, base='mean', limit=2):
	coord, data = BamCov(incov).to_matrix()
	if base=='mean':
		total = data.mean(axis=0)
	elif base=='median':
		total = np.median(data, axis=0)
	print(total, file=sys.stderr)
	data = data / total
	f = open(incov+'.snps', 'w')
	for (chrom, start, end), covs in zip(coord, data):
		if fold != 1:
			covs = covs * fold
#		covs = [min(round(cov,2), 2) for cov in covs]
#		covs = map(reset_cov, covs)
		covs = [reset_cov(cov, limit=limit) for cov in covs]

		id = '{}:{}'.format(chrom, start)
		line = [id, '0', '1'] + list(covs)
		line = map(str, line)
		line = '\t'.join(line)
		print(line, file=outgeno)
		#chrom = chrom.replace('chr', '')
		line = [id, start+1, chrom]
		line = map(str, line)
		line = ', '.join(line)
		print(line, file=f)
	f.close()
def reset_cov(cov, limit=2):
	return min(round(cov,2), limit)
def cov2gemma0(incov, outgeno):
	coord, data = BamCov(incov).to_matrix()
	for (chrom, start, end), covs in zip(coord, data):
		_len = end-start
		covs = covs / _len
		covs = map(define_gt, covs)
		id = '{}:{}'.format(chrom, start)
		line = [id, '0', '1'] + list(covs)
		line = map(str, line)
		line = '\t'.join(line)
		print(line, file=outgeno)
	
def define_gt(val):
	if val < 0.5:
		return 0
	elif val < 3:
		return 1
	else:
		return 2
def compare_gt(gemout, pheno, incov, outgt=sys.stdout, n=1, ):
	pheno = BimBam(pheno).get_column(n)
	groups = [g for g in sorted(set(pheno)) if g !='NA']
	#print(pheno, len(pheno), groups, file=sys.stderr)
	g0_idx = [i for i, p in enumerate(pheno) if p==groups[0]]
	g1_idx = [i for i, p in enumerate(pheno) if p==groups[1]]
#	sig_kmers = {rc.rs for rc in Gemma(gemout)}
	sig_kmers = {rc.split()[1] for rc in open(gemout)}
	print(len(sig_kmers), 'sig kmers', file=sys.stderr)
	opxs = []
	i = 0
	for line in BamCov(incov, ):
		i += 1
		if i % 100000 == 0:
			print(i, 'processed', file=sys.stderr)
		if line.chrom not in sig_kmers:
			continue
		line = np.array(line.covs)
		g0_dat = line[g0_idx]
		g1_dat = line[g1_idx]
		g0_mean = np.mean(g0_dat)
		g1_mean = np.mean(g1_dat)
		if g0_mean > g1_mean:
			op = '0>1'
		else:
			op = '0<1'
		opxs += [op]
	print( Counter(opxs))

def addgt(gemout, pheno, incov, outgt=sys.stdout, n=1, ):
	pheno = BimBam(pheno).get_column(n)
	groups = [g for g in sorted(set(pheno)) if g !='NA']
	g0_idx = [i for i, p in enumerate(pheno) if p==groups[0]]
	g1_idx = [i for i, p in enumerate(pheno) if p==groups[1]]
	coord, data = BamCov(incov, ).to_matrix()
	#data
	#print(data)
	#print(len(g0_idx), len(g1_idx))
	g0_dat = data[:, g0_idx]
	g1_dat = data[:, g1_idx]
	d_gt = {}
	for (chrom, start, end), g0, g1 in zip(coord, g0_dat, g1_dat):
		id = chrom
		d_gt[id] = (g0, g1)
	
	i = 0
	for rc in Gemma(gemout):		
		i += 1
		if i == 1:
			title = [rc.title] + ['sum0', 'sum1'] + ['0_{}'.format(j+1) for j in g0_idx] + ['1_{}'.format(j+1) for j in g1_idx]
			line = '\t'.join(title)
			print(line, file=outgt)
		if rc.rs not in d_gt:
			continue
		g0, g1 = d_gt[rc.rs]
		sum0, sum1 = sum(g0), sum(g1)
		line = [rc.line] + [sum0, sum1]+ list(g0) + list(g1)
		line = map(str, line)
		line = '\t'.join(line)
		print(line, file=outgt)

def cov2fold(incov, pheno, outpre=None, n=1, base='mean', win_size=500, win_step=None, 
		figfmt='pdf', normalize=1, chrs='', hw_ratio=None, marker=None, scatter=0):
	if win_step is None:
		win_step = win_size//2
	if outpre is None:
		outpre = '{}.{}'.format(pheno, n)
	chroms = chrs.split()
	pheno = BimBam(pheno).get_column(n)
	groups = [g for g in sorted(set(pheno)) if g !='NA']
	print(groups)
	g0_idx = [i for i, p in enumerate(pheno) if p==groups[0]]
	g1_idx = [i for i, p in enumerate(pheno) if p==groups[1]]
	print(g0_idx, g1_idx)
	coord, data = BamCov(incov, chroms).to_matrix()
	if normalize:
		if base=='mean':
			_mean = data.mean(axis=0)
		elif base=='median':
			_mean = np.median(data, axis=0)
	data = data / _mean
	g0_dat = data[:, g0_idx]
	g1_dat = data[:, g1_idx]
	g0_mean = g0_dat.mean(axis=1)
	g1_mean = g1_dat.mean(axis=1)
	fold = g1_mean / g0_mean

	fout = open(outpre + '.data', 'w')
	line = ['chrom', 'start', 'end', 'depth_{0}'.format(*groups), 'depth_{1}'.format(*groups), 'fold_{1}/{0}'.format(*groups)]
	print('\t'.join(line), file=fout)
	offset, last_offset, last_chr, last_start = 0, 0, 0, 0
	Xs, offsets, labels = [], [], []
	d_max = {}
	for (chrom, start, end), g0, g1, _fold in zip(coord, g0_mean, g1_mean, fold):
		if last_chr and last_chr != chrom:
			offset += last_start
			offsets += [offset]	# vlines
			labels += [last_chr]	# x labels
			d_max[last_chr] = last_start
		Xs += [start+ offset]
		last_chr = chrom
		last_start = start
		line = [chrom, start, end, g0, g1, _fold]
		line = map(str, line)
		print('\t'.join(line), file=fout)
	fout.close()
	offset += last_start
	offsets += [offset]
	labels += [last_chr]	# x labels
	d_max[last_chr] = last_start
	outfig = outpre + '.' + figfmt
	if hw_ratio is None:
		hw_ratio = sum(d_max.values()) / max(d_max.values())
	Xs, g0_mean, g1_mean, fold = bin_data(Xs, g0_mean, g1_mean, fold, win_size=win_size, win_step=win_step)
	fold = np.array(g1_mean) / np.array(g0_mean)
	bin_plot([Xs, Xs], outfig, Ys=[[g0_mean, g1_mean], [fold]], xlabels=labels, vlines=offsets, ylims=[3, 3], 
			ylabels=['mean depth', 'ratio'], hlines=[[0.5,1,2],[0.5,1,2]], hw_ratio=hw_ratio, marker=marker, scatter=scatter)
	outfig = outpre + '.ind.' + figfmt
	g0_dat = list(bin_data(*list(g0_dat.transpose()), win_size=win_size, win_step=win_step))
	g1_dat = list(bin_data(*list(g1_dat.transpose()), win_size=win_size, win_step=win_step))
	#print(g0_dat)
	bin_plot([Xs, Xs], outfig, Ys=[g0_dat, g1_dat], xlabels=labels, vlines=offsets, ylims=[3, 3],
			ylabels=['mean depth']*2, hlines=[[0.5,1,2],[0.5,1,2]], hw_ratio=hw_ratio)

def bin_data(*array, **args):
	for arr in array:
		yield _bin_data(arr, **args)
def _bin_data(arr, win_size=50, win_step=25):
	if win_size <=1:
		return  arr
	else:
		binned = []
		for i in range(0, len(arr), win_step):
			data = arr[i:i+win_size]
			binned += [np.mean(data)]
		return binned
def vlines2pos(vlines):
	bpos = []
	last_start = 0
	for lpos in vlines:
		bpos += [(last_start+lpos)/2]
		last_start = lpos
	return bpos
def bin_plot(Xs, outplot, Ys=[], xlabels=[], xlabel_positions=None, ylabels=[], vlines=[], vlines2=[], ylims=[], 
		hlines=[], hw_ratio=5, xlim=None, alpha=1, ls=None, marker=None, scatter=False):
	import matplotlib.pyplot as plt
	import matplotlib
	matplotlib.rcParams['xtick.minor.visible'] = True
	matplotlib.rcParams['ytick.minor.visible'] = True

	nsubplot = len(Ys)
	if not ylabels:
		ylabels = [None] * nsubplot
	if not ylims:
		ylims = [None] * nsubplot
	if not hlines:
		hlines = [None] * nsubplot
	if ls is None or isinstance(ls, str):
		ls = [ls] * nsubplot
	if marker is None or isinstance(marker, str):
		marker = [marker] * nsubplot
	if xlim is None:
		xlim = max(vlines)
	if xlabel_positions is None:
		xlabel_positions = vlines2pos(vlines)
	plt.figure(figsize=(8*hw_ratio, 5*nsubplot))
	i = 0
	for _Xs, _Ys, ylabel, ylim, hline, _ls, _marker in zip(Xs, Ys, ylabels, ylims, hlines, ls, marker):
		i +=1
		plt.subplot(nsubplot,1,i)
		for Y in _Ys:
			if scatter:
				plt.scatter(_Xs, Y, alpha=alpha, marker=_marker)
			else:
				plt.plot(_Xs, Y, alpha=alpha, ls=_ls, )
			if ylim is None:
				ylim = max(Y)
			plt.ylim(0, ylim)
			for v in vlines:
				plt.axvline(v, color="grey", )
			for v in vlines2:
				plt.axvline(v, color="red", )
		if ylabel is not None:
			plt.ylabel(ylabel)
		if hline is not None:
			if isinstance(hline, (int, float)):
				plt.axhline(hline, color="grey", ls='--')
			else:
				for h in hline:
					plt.axhline(h, color="grey", ls='--')
		plt.xlim(0, xlim)
	for xpos, label in zip(xlabel_positions, xlabels):
		y = -ylim/15
		plt.text(xpos, y, label, horizontalalignment='center',verticalalignment='top',fontsize=15)
	plt.subplots_adjust(hspace=0.02)
	plt.savefig(outplot, bbox_inches='tight')

def main():
	subcmd = sys.argv[1]
	kargs = parse_key_opts(sys.argv)
	if subcmd.startswith('cov2gemma'):
		incov = sys.argv[2]
		outgeno = sys.stdout
		if subcmd == 'cov2gemma':
			cov2gemma(incov, outgeno, **kargs)
		elif subcmd == 'cov2gemma0':
			cov2gemma0(incov, outgeno, **kargs)
	elif subcmd == 'cov2fold':
		incov = sys.argv[2]
		pheno = sys.argv[3]
		cov2fold(incov, pheno, **kargs)
	elif subcmd == 'addgt':
		gemout = sys.argv[2]
		pheno = sys.argv[3]
		incov = sys.argv[4]
		addgt(gemout, pheno, incov, **kargs)
	elif subcmd == 'cmpgt':
		gemout = sys.argv[2] 
		pheno = sys.argv[3]
		incov = sys.argv[4]
		compare_gt(gemout, pheno, incov, **kargs)
	elif subcmd == 'order_cov':
		incov, old_sd, new_sd = sys.argv[2:5]
		outcov = sys.stdout
		reorder_cov(incov, old_sd, new_sd, outcov)
	else:
		raise ValueError('Unknown subcmd: {}'.format(subcmd))

def parse_key_opts(args):
	d = {}
	pops = []
	for i, arg in enumerate(args):
		kv = arg.split('=', 1)
		if len(kv) != 2:
			continue
		pops += [i]
		key, val = kv
		val = tr_numeric(val)
		d[key] = val
	for i in sorted(pops, reverse=1):
		args.pop(i)
	return d

if __name__ == '__main__':
	main()
