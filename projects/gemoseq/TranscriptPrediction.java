package projects.gemoseq;

import java.io.File;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.TreeSet;
import java.util.stream.Collectors;

import de.jstacs.DataType;
import de.jstacs.parameters.EnumParameter;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterException;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.utils.IntList;
import projects.gemoma.Analyzer;
import projects.gemoma.ExtractRNAseqEvidence.Stranded;
import projects.gemoma.GeMoMaAnnotationFilter;
import projects.gemoseq.SplicingGraph.Gene;
import projects.gemoseq.SplicingGraph.Transcript;
import projects.gemoma.ReadStats;

public class TranscriptPrediction implements JstacsTool {
	
	
	private static class TranscriptResult{
		
		private LinkedList<Transcript> list;
		private SplicingGraph sg;
		
		public TranscriptResult(LinkedList<Transcript> list, SplicingGraph sg) {
			this.list = list;
			this.sg = sg;
		}
		
	}
	
	private static class Result implements Comparable<Result> {
		
		private int idx;
		private LinkedList<TranscriptResult> list;
		
		public Result(int idx, LinkedList<TranscriptResult> list) {
			this.idx = idx;
			this.list = list;
		}

		@Override
		public int compareTo(Result o) {
			return Integer.compare(this.idx, ((Result)o).idx);
		}
		
		
	}
	
	private static class OutputSet{
		
		private TreeSet<Result> internal;
		private int lastIdx;
		private String lastChrom;
		private int n;
		
		
		public OutputSet() {
			internal = new TreeSet<Result>();
			lastIdx = -1;
			lastChrom = "";
			this.n = 1;
		}
		
		
		public synchronized void add(int idx, LinkedList<TranscriptResult> list) {
			this.internal.add(new Result(idx,list));
		}
		
		public synchronized void print(BAMReader reader, PrintWriter wr, double minReadsPerGene, int minProteinLength, String genePrefix, boolean useChrPrefix) {
			
			while( this.internal.size() > 0 && this.internal.first().idx == lastIdx + 1) {

				Result res = this.internal.pollFirst();

				LinkedList<TranscriptResult> topList = res.list;

				for(TranscriptResult tres : topList) {

					SplicingGraph sg = tres.sg;
					LinkedList<Gene> genes = sg.finalize(tres.list, genePrefix,useChrPrefix, n, minReadsPerGene,minProteinLength);
					if(genes.size() > 0) {
						String chrom = genes.get(0).getChrom();
						if(!lastChrom.contentEquals(chrom)) {
							int len = reader.getSequenceLength(chrom);
							wr.println("##sequence-region "+chrom+" 1 "+len);
						}


						lastChrom = chrom;

					}

					n += genes.size();
					for(Gene g : genes) {
						wr.print(g);
					}

					wr.flush();

				}
				
				lastIdx = res.idx;

			}
		}
		
	}
	
	private static class Config{
		private int minIntronLength;
		private int maxIntronLength;
		private Stranded stranded; 
		private double minReads; 
		private double minFraction;
		private double minIntronReads;
		private double minIntronFraction;
		private double maxNumTranscripts; 
		private double percentExplained; 
		private double minReadsPerTranscript;
		private double minReadsPerGene;
		private double maxFraction;
		private double percentAbundance;
		private double scaleIntronReads;
		private double delta;
		private int nIterations;
		private ReadStats stats;
		private int minProteinLength;
		private double spilloverFactor;
		private int maxGapFilled;
		private int maxMM;
		
		private boolean doSplits;
		private boolean filterSpillover;
		private boolean longReads;
		
		private String geneBase;
		private boolean useChrPrefix;
		
		public Config(int minIntronLength, int maxIntronLength, Stranded stranded, double minReads, double minFraction,
				double minIntronReads, double minIntronFraction, double maxNumTranscripts, double percentExplained,
				double minReadsPerTranscript, double minReadsPerGene, double maxFraction, double percentAbundance, double scaleIntronReads,
				double delta, int nIterations, ReadStats stats, int minProteinLength, int maxGapFilled, boolean longReads, String geneBase, boolean useChrPrefix) {
			this.minIntronLength = minIntronLength;
			this.maxIntronLength = maxIntronLength;
			this.stranded = stranded;
			this.minReads = minReads;
			this.minFraction = minFraction;
			this.minIntronReads = minIntronReads;
			this.minIntronFraction = minIntronFraction;
			this.maxNumTranscripts = maxNumTranscripts;
			this.percentExplained = percentExplained;
			this.minReadsPerTranscript = minReadsPerTranscript;
			this.minReadsPerGene = minReadsPerGene;
			this.maxFraction = maxFraction;
			this.percentAbundance = percentAbundance;
			this.scaleIntronReads = scaleIntronReads;
			this.delta = delta;
			this.nIterations = nIterations;
			this.stats = stats;
			this.minProteinLength = minProteinLength;
			this.spilloverFactor = 10;
			this.maxGapFilled = maxGapFilled;
			this.longReads = longReads;
			if(longReads) {
				this.maxMM = Integer.MAX_VALUE;
			}else {
				this.maxMM = 3;
			}
			
			this.doSplits = true;
			this.filterSpillover = true;
			this.geneBase = geneBase;
			this.useChrPrefix = useChrPrefix;
			
		}
		
		
		
	}
	
	private enum State{
		IDLE,
		COMP,
		PRINT,
		EXCEPTION
	}
	
	
	private static class WorkerPool{
		
		private Worker[] workers;
		private Config config;
		
		public WorkerPool(Config config, int numThreads) {
			this.config = config;
			workers = new Worker[numThreads];
			for(int i=0;i<workers.length;i++) {
				workers[i] = new Worker(this,config,i);
				(new Thread(workers[i])).start();
			}
		}
		
		public void stopAll() {
			while(true) {
				int stopped = 0;
				for(int i=0;i<workers.length;i++) {
					if(workers[i].getState() == State.IDLE || workers[i].getState() == State.EXCEPTION) {
						workers[i].run = false;
						synchronized(workers[i]) {
							workers[i].notify();
						}
						stopped++;
					}
				}
				
				if(stopped == workers.length) {
					return;
				}
				
				try {
					synchronized(this) {
						wait(1000);
					}
				} catch (InterruptedException e) {
				}
				
			}
		}
		
		public synchronized Worker getFree() throws Exception {
			
			while(true) {
				
				for(int i=0;i<workers.length;i++) {
					if(workers[i].getState() == State.EXCEPTION) {
						throw workers[i].exception;
					}
				}
				
				for(int i=0;i<workers.length;i++) {

					if(workers[i].getState() == State.IDLE) {
						return workers[i];
					}
				}
				
				try {
					wait(1000);
				} catch (InterruptedException e) {
				}
				
			}
			
		}
		
	}
	
	private static class Worker implements Runnable {

		private WorkerPool pool;
		private boolean run;
		private Config config;
		private State state;
		private BAMReader reader;
		private PrintWriter wr;
		private OutputSet outset;
		private int idx;
		private Region region;
		private Exception exception;
		
		private int id;
		
		public Worker(WorkerPool pool, Config config, int id) {
			this.pool = pool;
			this.run = true;
			this.config = config;
			this.state = State.IDLE;
			this.id = id;
		}
		
		
		public synchronized State getState() {
			return state;
		}


		public synchronized void setPrint(BAMReader reader, PrintWriter wr, OutputSet outset) {
			this.reader = reader;
			this.wr = wr;
			this.outset = outset;
			this.state = State.PRINT;
			this.notify();
		}
		
		public synchronized void setCompute(int idx, Region region, OutputSet outset) {
			this.idx = idx;
			this.region = region;
			this.outset = outset;
			this.state = State.COMP;
			this.notify();
		}
		
		
		@Override
		public void run() {
			while(run) {
				
				try {
				synchronized(this) {
					if(state == State.IDLE) {
						try {
							wait(1000);
						} catch (InterruptedException e) {
							e.printStackTrace();
							return;
						}
					}
				}
				
				if(state == State.COMP) {
					
					try {
						compute(idx,region,outset);
					} catch (CloneNotSupportedException e) {
						e.printStackTrace();
						return;
					}
					
					state = State.IDLE;
					
					synchronized(pool) {
//						System.out.println("notify");
//						System.out.flush();
						pool.notify();
					}
					
				}else if(state == State.PRINT) {
					
					outset.print(reader, wr, config.minReadsPerGene,config.minProteinLength,config.geneBase, config.useChrPrefix);
					
					state = State.IDLE;
					
					synchronized(pool) {;
						pool.notify();
					}
				}
				
				}catch(Exception e) {
					e.printStackTrace();
					exception = e;
					state = State.EXCEPTION;
					synchronized(pool) {;
						pool.notify();
					}
					return;
				}
				
				
			}
		}
		
		private void compute(int idx, Region region, OutputSet outset) throws CloneNotSupportedException {
			
			boolean filterBefore = true;
			
			LinkedList<TranscriptResult> topList = new LinkedList<TranscriptResult>();
			
			double rpk = region.getTheoreticalNumberOfReads()/(double)(region.getRegionEnd()-region.getRegionStart()+1)*1000.0;

			ReadGraph rg2 = region.buildGraph(config.minIntronLength, config.maxGapFilled, config.stranded, config.maxMM, config.longReads);			
			
			rg2.pruneByAbsoluteNumberOfReads(config.minReads,false);
			if(rpk > 200) {
				rg2.pruneByAbsoluteNumberOfReads(config.minIntronReads*2,true);
			}else {
				rg2.pruneByAbsoluteNumberOfReads(config.minIntronReads,true);
			}
			
			rg2.pruneByRelativeNumberOfReads(config.minIntronFraction,true);
			
			
			rg2.pruneBySplitLength(config.stats);
			
			if(config.longReads) {
				rg2.pruneAlternativeIntronsLong(5);
			}
			
			ReadGraph rg = null;
			LinkedList<ReadGraph.Edge> intronEdges = rg2.getIntronEdges().stream().map(e -> new ReadGraph.Edge(e.getRelStart(),e.getRelEnd(),e.getNumberOfReads()) ).collect(Collectors.toCollection(LinkedList::new));
			int edgeRegionStart = rg2.regionStart;
			while( (rg = rg2.nextSplit()) != null ) {

				rg.condense();

				int[] proposedSplits = new int[0];
				
				if(config.doSplits) {
					proposedSplits = rg.proposeSplitByReadCoverage(200, 20);
				}

				SplicingGraph sg = new SplicingGraph(rg,region,config.minIntronLength,proposedSplits,config.longReads);
				
				LinkedList<Transcript> list = new LinkedList<SplicingGraph.Transcript>();
				if(config.longReads) {
					sg.enumerateLongReadTranscripts(list,config.minReads);
				}else {
					sg.enumerateTranscripts2(list,config.minReads,config.minFraction,config.maxNumTranscripts);//TODO replace for long
				}
				
				list.stream().forEach(e -> {if(e.getStrand() == '.') e.setStrand(region.getStrand());} );
								
				int sizeBefore = list.size();
				
				LinkedList<Transcript> sList = list;
				
				
				if(config.doSplits && proposedSplits.length>0) {
					sList = sg.testSplitByORF(sList,proposedSplits,config.minProteinLength);
					
				}
				
				if(config.doSplits && config.stranded == Stranded.FR_UNSTRANDED) {
					sList = sg.testSplitBySpliceSitesAndORF(sList,config.minProteinLength);//TODO before first heuristic?
					
				}
				
				
				if(config.doSplits && sList.size() == sizeBefore) {
					sList = sg.testSplitByORF(sList,config.minProteinLength);
				}
				
				LinkedList<Transcript>[] sList2 = sg.splitByLocationAndMakeUnique(sList,config.minProteinLength,config.stranded);
				
				
				for(int i=0;i<sList2.length;i++) {
				
					
					if(config.filterSpillover && filterBefore) {
						sList2[i] = filterSpilloverTranscripts(sList2[i],region,config.spilloverFactor);
						sList2[i] = filterSpilloverIntronTranscripts(sList2[i],intronEdges, edgeRegionStart, config.spilloverFactor);
					}
					
					list = sg.quantify(sList2[i], config.percentExplained, config.minReadsPerTranscript, config.maxFraction, config.percentAbundance, 
							config.scaleIntronReads, config.delta, config.nIterations, region.getScale(), config.minIntronLength, config.longReads);
					
					list = filterImplausibleTranscripts(list,config.stranded);
					list = filterSingleExonGenes(list,config.minReadsPerTranscript);
					if(config.filterSpillover && !filterBefore) {
						list = filterSpilloverTranscripts(list,region,config.spilloverFactor);
						list = filterSpilloverIntronTranscripts(list,intronEdges, rg2.regionStart, config.spilloverFactor);						
					}
					topList.add(new TranscriptResult(list, sg));
				}			
			}
			outset.add(idx, topList);
		}

		
		private LinkedList<Transcript> filterSingleExonGenes(LinkedList<Transcript> list, double minReadsPerTranscript) {
			Iterator<Transcript> it = list.iterator();
			IntList li = new IntList();
			int i=0;
			while(it.hasNext()) {
				Transcript t = it.next();
				if(t.getNumberOfIntrons()==0) {
					
					double abundance = t.getAbundance();
					
					
					if(abundance*(1.0-t.getUTRFraction()) < minReadsPerTranscript) {
						li.add(i);
					}
					
				}
				i++;
			}
			for(int j=li.length()-1;j>=0;j--) {
				list.remove( li.get(j) );
			}

			return list;
		}


		private LinkedList<Transcript> filterSpilloverIntronTranscripts(LinkedList<Transcript> list, LinkedList<ReadGraph.Edge> intronEdges, int edgeRegionStart, double factor) {
			Iterator<Transcript> it = list.iterator();
			IntList li = new IntList();
			int i=0;
			while(it.hasNext()) {
				Transcript t = it.next();
				if(t.getNumberOfIntrons()==0 || t.percentCanonicalSpliceSites()<0.75) {
					
					double myCov = t.getMaximumNumberOfReads();
					double intronCov = t.getMaximumNumberOfCoveringIntronReads(intronEdges, edgeRegionStart);
					
					if(myCov*factor*(1.0-t.getUTRFraction()) < intronCov) {
						li.add(i);
					}
					
				}
				i++;
			}
			for(int j=li.length()-1;j>=0;j--) {
				list.remove( li.get(j) );
			}

			return list;
		}
		

		private LinkedList<Transcript> filterSpilloverTranscripts(LinkedList<Transcript> list, Region region2, double factor) {
			if(region2.getRevRegion() != null) {
				Iterator<Transcript> it = list.iterator();
				IntList li = new IntList();
				int i=0;
				while(it.hasNext()) {
					Transcript t = it.next();
					if(t.getNumberOfIntrons()<1 || t.percentCanonicalSpliceSites() < 0.75) {//TODO?

						Region revRegion = region2.getRevRegion();
						
						int[] covRegion = region2.getCoverageByBlocks();
						int[] covRev = revRegion.getCoverageByBlocks();

						double sumRegion = 0;
						double sumRev = 0;
						for(int j=Math.max(t.getStart(),region2.getRegionStart());j<t.getEnd() && j-region2.getRegionStart()<covRegion.length;j++) {
							sumRegion += covRegion[ j-region2.getRegionStart() ];
							if(j-revRegion.getRegionStart() >= 0 && j-revRegion.getRegionStart()<covRev.length) {
								sumRev += covRev[ j-revRegion.getRegionStart() ];
							}
						}
						if(sumRegion*factor*(1.0-t.getUTRFraction())<sumRev) {
							li.add(i);
						}
					}
					i++;
				}
				for(int j=li.length()-1;j>=0;j--) {
					list.remove( li.get(j) );
				}
			}
			return list;
		}


		private LinkedList<Transcript> filterImplausibleTranscripts(LinkedList<Transcript> linkedList,
				Stranded stranded) {
			if(stranded != Stranded.FR_UNSTRANDED) {
				Iterator<Transcript> it = linkedList.iterator();
				IntList li = new IntList();
				int i=0;
				while(it.hasNext()) {
					Transcript t = it.next();
					char strand = t.getStrand();
					if(strand != '.') {
						
						char intronStrand = t.getStrandByIntrons();
						if(intronStrand != '.' && intronStrand != strand) {
							li.add(i);
						}
					}
					i++;
				}
				for(int j=li.length()-1;j>=0;j--) {
					linkedList.remove( li.get(j) );
				}
			}
			return linkedList;
		}
		
	}
	
	
	

	public static void main(String[] args) throws Exception {
		
		CLI cli = new CLI(new boolean[] {true,false,false,false,false},new TranscriptPrediction(), new PredictCDSFromGFF(), new GeMoMaAnnotationFilter(), new Analyzer(), new MergeGeMoMaGeMoSeq());
		
		cli.run(args);
	}
	
	
	public TranscriptPrediction() {
	}

	@Override
	public ToolParameterSet getToolParameters() {
		
		LinkedList<Parameter> pars = new LinkedList<Parameter>();
		
		pars.add(new FileParameter("Genome", "Genome sequence as FastA", "fa,fna,fasta", true));
		pars.add(new FileParameter("Mapped reads","Mapped Reads in BAM format, coordinate sorted","bam",true));
		
		try {
			pars.add(new EnumParameter(Stranded.class, "Library strandedness", true,Stranded.FR_UNSTRANDED.name()));
			
			
			pars.add(new SimpleParameter(DataType.INT,"Longest intron length","Length of the longest intron reported",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),100000));
			pars.add(new SimpleParameter(DataType.INT,"Shortest intron length","Length of the shortest intron considered",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),10));
			
			pars.add(new SimpleParameter(DataType.BOOLEAN,"Long reads","Long-read mode",true,false));
		
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum number of reads","Minimum number of reads required for an edge in the read graph",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),1.0));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum fraction of reads","Minimum fraction of reads relative to adjacent exons that must support an intron in the enumeration",true,new NumberValidator<Double>(0.0, 1.0),0.01));
			
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum number of intron reads","Minimum number of reads required for an intron",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),1.0));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum fraction of intron reads","Minimum fraction of reads relative to adjacent exons for an intron to be considered",true,new NumberValidator<Double>(0.0, 1.0),0.01));
			
			pars.add(new SimpleParameter(DataType.DOUBLE,"Percent explained","Percent of abundance that must be explained by transcript models after quantification",true,new NumberValidator<Double>(0.0, 1.0),0.9));//was 0.8
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum reads per gene","Minimum abundance required for a gene to be reported",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),40.0));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Minimum reads per transcript","Minimum abundance required for a transcript to be reported",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),20.0));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Percent abundance","Minimum relative abundance required for a transcript to be reported",true,new NumberValidator<Double>(0.0, 1.0),0.05));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Successive fraction","Factor of the drop in abundance between successive transcript models",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),20.0));//was 10.0
			
			pars.add(new SimpleParameter(DataType.INT,"Maximum region length","Maximum length of a region considered before it is split",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),750000));
			pars.add(new SimpleParameter(DataType.DOUBLE,"Maximum region coverage","Maximum coverage in a region before reads are down-sampled",true,new NumberValidator<Double>(0.0, Double.POSITIVE_INFINITY),100.0));
			pars.add(new SimpleParameter(DataType.INT,"Maximum filled gap length","Maximum length of a gap filled by dummy reads",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),50));
			pars.add(new SimpleParameter(DataType.INT,"Quality filter","Minimum mapping quality required for a read to be considered",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),40));
			
			pars.add(new SimpleParameter(DataType.INT,"Minimum protein length","Minimum length of protein in AA",true,new NumberValidator<Integer>(0, Integer.MAX_VALUE),70));
			
			pars.add(new SimpleParameter(DataType.STRING, "Gene prefix", "Prefix to add to all gene names", true,"G"));
			pars.add(new SimpleParameter(DataType.BOOLEAN, "Gene names with chromosome", "If true, gene names will be constructed as <Gene prefix><chr>.<geneNumber>. Gene numbers will be assigned successively across all chromosomes.", true, false));
						
		}catch(ParameterException ex) {
			ex.printStackTrace();
		}
		
		return new ToolParameterSet(this.getToolName(),pars);
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads)
			throws Exception {
		
		WorkerPool pool = null;
		
		try {
			if(threads>6) {
				protocol.appendWarning("LARGE NUMBER OF THREADS DETECTED!\nCurrently, I/O of GeMoSeq runs on a single thread and runtime is limited by I/O performance. Hence, running GeMoSeq with a large number of threads is not recommended. On our infrastructure, a number of 6 threads has been the sweet spot.\n\n");
			}

			String genome = ((FileParameter)parameters.getParameterForName("Genome")).getFileContents().getFilename();
			String bamFile = ((FileParameter)parameters.getParameterForName("Mapped reads")).getFileContents().getFilename();;
			int minIntronLength = (int) parameters.getParameterForName("Shortest intron length").getValue();
			int maxIntronLength = (int) parameters.getParameterForName("Longest intron length").getValue();

			double maxCov = (double) parameters.getParameterForName("Maximum region coverage").getValue();//100.0;
			double sample = 0.5;


			Stranded stranded = (Stranded) ((EnumParameter)parameters.getParameterForName(Stranded.class.getSimpleName())).getValue();
			double minReads = (double) parameters.getParameterForName("Minimum number of reads").getValue();
			double minFraction = (double) parameters.getParameterForName("Minimum fraction of reads").getValue();
			double minIntronReads = (double) parameters.getParameterForName("Minimum number of intron reads").getValue();
			double minIntronFraction = (double) parameters.getParameterForName("Minimum fraction of intron reads").getValue();

			double maxNumTranscripts = 1000.0;
			double percentExplained = (double) parameters.getParameterForName("Percent explained").getValue();
			double minReadsPerTranscript = (double) parameters.getParameterForName("Minimum reads per transcript").getValue();
			double maxFraction = (double) parameters.getParameterForName("Successive fraction").getValue();
			double percentAbundance = (double) parameters.getParameterForName("Percent abundance").getValue();
			double scaleIntronReads = 5.0;//TODO was 2.0
			double delta = 1E-4;
			int nIterations = 10000;
			double minReadsPerGene = (double) parameters.getParameterForName("Minimum reads per gene").getValue();

			int minQuality = (int) parameters.getParameterForName("Quality filter").getValue();
			int maxLen = (int) parameters.getParameterForName("Maximum region length").getValue();
			int maxGap = (int) parameters.getParameterForName("Maximum filled gap length").getValue();

			int minProteinLength = (int) parameters.getParameterForName("Minimum protein length").getValue();

			boolean longReads = (boolean) parameters.getParameterForName("Long reads").getValue();

			String geneBase = (String) parameters.getParameterForName("Gene prefix").getValue();
			boolean useChrPrefix = (boolean) parameters.getParameterForName("Gene names with chromosome").getValue();

			ReadStats stats = new ReadStats(minIntronLength, 1.0, bamFile);


			Config config = new Config(minIntronLength, maxIntronLength, stranded, minReads, minFraction, minIntronReads, minIntronFraction, 
					maxNumTranscripts, percentExplained, minReadsPerTranscript, minReadsPerGene, maxFraction, percentAbundance, scaleIntronReads, 
					delta, nIterations,stats,minProteinLength,maxGap,longReads,geneBase, useChrPrefix);


			BAMReader reader = new BAMReader(maxIntronLength, bamFile, maxCov, sample, stranded,minQuality,maxLen, maxGap, longReads);





			Genome.init(genome);


			File out = File.createTempFile("predictions", ".temp", new File("."));
			out.deleteOnExit();

			PrintWriter wr = new PrintWriter(out);

			OutputSet outset = new OutputSet();

			
			Worker worker = null;
			if(threads > 1) {
				pool = new WorkerPool(config, threads-1);

			}else {
				worker = new Worker(null,config,0);
			}

			int idx = 0;

			wr.println("##gff-version 3");
			while(reader.hasNext()) {

				Region region2 = reader.next();

				Region[] regions = new Region[] {region2};


				if(region2.getRegionEnd()-region2.getRegionStart()>maxLen) {
					regions = region2.splitByCov(maxLen);
					System.out.println("#split: "+region2.getChrom()+" "+Arrays.toString(regions));
					System.out.flush();
				}

				for(Region region : regions) {

					if(pool == null) {
						worker.compute(idx, region, outset);
						outset.print(reader, wr, minReadsPerGene,config.minProteinLength,config.geneBase,config.useChrPrefix);
					}else {

						worker = pool.getFree();
						worker.setCompute(idx, region, outset);

						worker = pool.getFree();
						worker.setPrint(reader, wr, outset);

					}

					idx++;
				}

			}

			if(pool != null) {
				pool.stopAll();
			}

			outset.print(reader, wr, minReadsPerGene,config.minProteinLength,config.geneBase,config.useChrPrefix);

			wr.close();

			TextResult tr = new TextResult("Transcript Predictions", "Predictions in GFF format", new FileParameter.FileRepresentation(out.getAbsolutePath()), "gff3", getToolName(), null, true);

			return new ToolResult("Result of "+getToolName(), getToolName(), null, new ResultSet(tr), parameters, getToolName(), new Date(System.currentTimeMillis()) );
		}catch(Throwable e) {
			if(pool != null) {
				pool.stopAll();
			}
			throw e;
		}
	}

	@Override
	public String getToolName() {
		return "GeMoSeq";
	}

	@Override
	public String getToolVersion() {
		return "1.2.3";
	}

	@Override
	public String getShortName() {
		return "gemoseq";
	}

	@Override
	public String getDescription() {
		return "transcript reconstruction from RNA-seq data";
	}

	@Override
	public String getHelpText() {
		return "";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		return null;
	}

	@Override
	public void clear() {

	}

	@Override
	public String[] getReferences() {
		return null;
	}

}
