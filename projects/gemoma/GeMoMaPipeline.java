package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintStream;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorCompletionService;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import de.jstacs.DataType;
import de.jstacs.io.FileManager;
import de.jstacs.io.RegExFilenameFilter;
import de.jstacs.io.RegExFilenameFilter.Directory;
import de.jstacs.io.XMLParser;
import de.jstacs.parameters.AbstractSelectionParameter;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.FileParameter.FileRepresentation;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.utils.IntList;
import de.jstacs.utils.SafeOutputStream;
import projects.FastaSplitter;

/**
 * The GeMoMa pipeline as one tool.
 * 
 * @author Jens Keilwagen
 * 
 * @see ExtractRNAseqEvidence
 * @see Extractor
 * @see GeMoMa
 * @see GAF
 */
public class GeMoMaPipeline implements JstacsTool {

	//XXX improve runtime? https://docs.oracle.com/javase/8/docs/api/java/util/concurrent/CountDownLatch.html

	static int maxSize;
	static long timeOut, maxTimeOut;
	
	static {
		try {
			File jarfile = new File(Galaxy.class.getProtectionDomain().getCodeSource().getLocation().toURI().getPath());
			File ini = new File( jarfile.getParentFile().getAbsolutePath() + File.separator + "GeMoMa.ini.xml" );
			StringBuffer xml = FileManager.readFile(ini);
			maxSize = XMLParser.extractObjectForTags(xml, "maxSize", Integer.class);
			timeOut = XMLParser.extractObjectForTags(xml, "timeOut", Long.class);
			maxTimeOut = XMLParser.extractObjectForTags(xml, "maxTimeOut", Long.class);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	String home;
	String target;
	int threads, gapOpen, gapExt;
	double eValue;
	boolean rnaSeq, stop;
	
	ToolParameterSet extractorParams, gemomaParams, gafParams;
	ExecutorService queue;
	ExecutorCompletionService ecs;
	
	ArrayList<Species> species;
	int speciesCounter = 0;
	IntList splitsPerSpecies;
	
	RNASeq rnaSeqData;
	
	ArrayList<Process> process;
	ArrayList<FlaggedRunnable> list;
	
	HashMap<String,String> selected;
	
	public ToolParameterSet getRelevantParameters( ToolParameterSet params, String... remove ) {
		HashSet<String> removeNames = new HashSet<String>();
		for( String r : remove ) {
			removeNames.add(r);
		}
		ArrayList<Parameter> list = new ArrayList<Parameter>();
		for( int i = 0; i < params.getNumberOfParameters(); i++ ) {
			Parameter p = params.getParameterAt(i);
			if( !removeNames.contains(p.getName()) ) {
				list.add(p);
			}			
		}
		return new ToolParameterSet( params.getToolName(), list.toArray(new Parameter[0])  );
	}
	
	private static String buscoPath= "BUSCO-references/";
	private static FilenameFilter dirFilter = new RegExFilenameFilter("", Directory.REQUIRED, true, ".*" );
	
	public ToolParameterSet getToolParameters() {
		ToolParameterSet ere = new ExtractRNAseqEvidence().getToolParameters();
		ToolParameterSet ex = getRelevantParameters(new Extractor(maxSize).getToolParameters(), "annotation", "genome", "selected", "verbose", "genetic code", Extractor.name[3], Extractor.name[4] );
		ToolParameterSet gem = getRelevantParameters(new GeMoMa(maxSize,timeOut,maxTimeOut).getToolParameters(), "tblastn results", "target genome", "cds parts", "assignment", "query proteins", "selected", "verbose", "genetic code", "tag", "coverage" );
		ToolParameterSet gaf = getRelevantParameters(new GeMoMaAnnotationFilter().getToolParameters(), "predicted annotation", "tag");
		
		ArrayList<String> keys = new ArrayList<String>();
		keys.add("own");
		keys.add("pre-extracted");
		
		try{
			ArrayList<SimpleParameterSet> values = new ArrayList<SimpleParameterSet>();
			values.add( new SimpleParameterSet(
							new SimpleParameter(DataType.STRING,"ID","ID to distinguish the different reference species", false, ""),
							new FileParameter( "annotation", "Reference annotation file (GFF or GTF), which contains gene models annotated in the reference genome", "gff,gtf", true ),
							new FileParameter( "genome", "Reference genome file (FASTA)", "fasta,fa.gz,fasta.gz",  true )
			) );
			values.add( new SimpleParameterSet(
							new SimpleParameter(DataType.STRING,"ID","ID to distinguish the different reference species", false, ""),
							new FileParameter( "cds parts", "The query cds parts file (FASTA), i.e., the cds parts that have been blasted", "fasta", true ),
							new FileParameter( "assignment", "The assignment file, which combines parts of the CDS to transcripts", "tabular", false ),
							new FileParameter( "query proteins", "optional query protein file (FASTA) for computing the optimal alignment score against complete protein prediction", "fasta", false )
			) );
			
			File busco = new File( buscoPath );
			if( busco.exists() && busco.isDirectory() ) {
				File[] d = busco.listFiles(dirFilter);
				
				String[] k = new String[d.length];
				for( int i = 0; i < d.length; i++ ) {
					k[i] = d[i].getName();
				}
				keys.add("busco");
				values.add( new SimpleParameterSet( new SelectionParameter(DataType.STRING, k, k, "BUSCO clades", "the clade to be used", true) ) );
			}
			
			return new ToolParameterSet( getShortName(),
				new FileParameter( "target genome", "Target genome file (FASTA)", "fasta,fa.gz,fasta.gz",  true ),

				new ParameterSetContainer( "reference species", "", 
						new ExpandableParameterSet( new SimpleParameterSet(	
								new SelectionParameter( DataType.PARAMETERSET, keys.toArray(new String[0]),
										values.toArray(new SimpleParameterSet[0]), "species", "data for reference species", true
								)
						), "reference", "", 1 )
				),/**/
				
				new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns can be used to determine a target region that should be overlapped by the prediction, if columns 2 to 5 contain chromosome, strand, start and end of region", "tabular,txt", maxSize>-1 ),
				new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
				new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, GeMoMa.TAG ),
				
				new SelectionParameter( DataType.PARAMETERSET, 
						new String[]{"NO", "MAPPED","EXTRACTED"},
						new Object[]{
								new SimpleParameterSet(),
								ere,
								new SimpleParameterSet(
										new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
												new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff", false )
											), "introns", "", 1 ) ),
	
										new ParameterSetContainer( "coverage", "", new ExpandableParameterSet( new SimpleParameterSet(	
											new SelectionParameter( DataType.PARAMETERSET, 
													new String[]{"NO", "UNSTRANDED", "STRANDED"},
													new Object[]{
														//no coverage
														new SimpleParameterSet(),
														//unstranded coverage
														new SimpleParameterSet(
																new FileParameter( "coverage_unstranded", "The coverage file contains the unstranded coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
														),
														//stranded coverage
														new SimpleParameterSet(
																new FileParameter( "coverage_forward", "The coverage file contains the forward coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true ),
																new FileParameter( "coverage_reverse", "The coverage file contains the reverse coverage of the genome per interval. Intervals with coverage 0 (zero) can be left out.", "bedgraph", true )
														)
													},  "coverage", "experimental coverage (RNA-seq)", true
											)
										), "coverage", "", 1 ) )
								)
						},
						"RNA-seq evidence", "data for RNA-seq evidence", true ),

				new ParameterSetContainer("Extractor parameter set", "parameters for the Extrator module of GeMoMa", ex ),
				new ParameterSetContainer("GeMoMa parameter set", "parameters for the GeMoMa", gem ),
				new ParameterSetContainer("GAF parameter set", "parameters for the GAF module of GeMoMa", gaf ),
				
				new SimpleParameter( DataType.BOOLEAN, "predicted proteins", "If *true*, returns the predicted proteins of the target organism as fastA file", true, true ),
				new SimpleParameter( DataType.BOOLEAN, "predicted CDSs", "If *true*, returns the predicted CDSs of the target organism as fastA file", true, false ),
				new SimpleParameter( DataType.BOOLEAN, "predicted genomic regions", "If *true*, returns the genomic regions of predicted gene models of the target organism as fastA file", true, false ),
				
				new SimpleParameter( DataType.BOOLEAN, "debug", "If *false* removes all temporary files even if the jobs exits unexpected", true, true )
			);		
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}
	
	/**
	 * This method is intended to set the parameters of GeMoMa modules as given in the global parameters.
	 * 
	 * @param given the given (global) {@link ParameterSet}
	 * @param toBeSet the {@link ParameterSet} of a GeMoMa module
	 * 
	 * @throws IllegalValueException should not happen
	 * @throws CloneNotSupportedException should not happen
	 */
	static void setParameters( ParameterSet given, ParameterSet toBeSet ) throws IllegalValueException, CloneNotSupportedException {
		if( given instanceof ExpandableParameterSet ) {
			int n = ((ExpandableParameterSet) given).getNumberOfParameters();
			ExpandableParameterSet eTBS = (ExpandableParameterSet) toBeSet;
			while( eTBS.getNumberOfParameters() < n ) {
				eTBS.addParameterToSet();
			}
		}
		for( int i = 0; i < given.getNumberOfParameters(); i++ ) {
			Parameter p = given.getParameterAt(i);
			Parameter q = toBeSet.getParameterForName( p.getName() );
			/*
			System.out.println(q.getName() + "\t" + p.getName() );
			System.out.println(q.getClass().getName() + "\t" + p.getClass().getName() );
			*/
			if( p instanceof ParameterSetContainer ) {
				setParameters( ((ParameterSetContainer)p).getValue(), ((ParameterSetContainer)q).getValue() );
			} else if( p instanceof FileParameter ){
				if( p.getValue() != null ) {
					q.setValue( ((FileParameter) p).getFileContents() );
				} else {
					q.reset();
				}
			} else if( p instanceof AbstractSelectionParameter && ((AbstractSelectionParameter)p).getDatatype()==DataType.PARAMETERSET ) {
				AbstractSelectionParameter a = (AbstractSelectionParameter) p;
				ParameterSet ps = (ParameterSet) a.getValue();
				AbstractSelectionParameter b = (AbstractSelectionParameter) q;
				//b.setValue(ps.getName());
				setParameters( ps, (ParameterSet) b.getValue() );
			} else {
				q.setValue( p.getValue() );
			}
		}
	}
	
	void setGeMoMaParams( Protocol protocol ) throws IllegalValueException, CloneNotSupportedException {
		gemomaParams.getParameterForName("target genome").setValue(target);

		//introns
		ExpandableParameterSet exp = (ExpandableParameterSet) ((ParameterSetContainer) gemomaParams.getParameterForName("introns")).getValue();
		for( int i = 0; i < rnaSeqData.introns.size(); i++ ) {
			if( exp.getNumberOfParameters() <= i ) {
				exp.addParameterToSet();
			}
			((ParameterSetContainer) exp.getParameterAt(i)).getValue().getParameterAt(0).setValue(rnaSeqData.introns.get(i));
		}

		//coverage
		exp = (ExpandableParameterSet) ((ParameterSetContainer) gemomaParams.getParameterForName("coverage")).getValue();
		int i = 0;
		for( int j = 0; j < rnaSeqData.coverageUn.size(); j++, i++ ) {
			if( exp.getNumberOfParameters() <= i ) {
				exp.addParameterToSet();
			}
			SelectionParameter sel = (SelectionParameter) ((SimpleParameterSet)(exp.getParameterAt(i).getValue())).getParameterAt(0);
			sel.setValue("UNSTRANDED");
			((SimpleParameterSet) sel.getValue()).getParameterAt(0).setValue(rnaSeqData.coverageUn.get(j));
		}
		for( int j = 0; j < rnaSeqData.coverageFwd.size(); j++, i++ ) {
			if( exp.getNumberOfParameters() <= i ) {
				exp.addParameterToSet();
			}
			SelectionParameter sel = (SelectionParameter) ((SimpleParameterSet)(exp.getParameterAt(i).getValue())).getParameterAt(0);
			sel.setValue("STRANDED");
			((SimpleParameterSet) sel.getValue()).getParameterAt(0).setValue(rnaSeqData.coverageFwd.get(j));
			((SimpleParameterSet) sel.getValue()).getParameterAt(2).setValue(rnaSeqData.coverageRC.get(j));
		}
		
		//CLI.print("", null, gemomaParams, "", protocol, "");;
	}

	ArrayList<Result> res;
	Protocol pipelineProtocol;
	
	int phase = 0;
	
	private void addNewPhase() {
		String s = "starting phase " + ++phase;
		int l = s.length();
		s+="\n";
		for( int i = 0; i < l; i++ ) {
			s += "=";
		}
		pipelineProtocol.append("\n" + s + "\n");
	}
	
	private class Species {
		String id, cds, assignment, protein = null;
		
		Species( String id ) {
			this.id = id;
		}
	}
	
	private class RNASeq {
		ArrayList<String> introns = new ArrayList<String>();
		ArrayList<String> coverageUn =  new ArrayList<String>();
		ArrayList<String> coverageFwd =  new ArrayList<String>();
		ArrayList<String> coverageRC =  new ArrayList<String>();
	}
	
	private static void setParameters( ToolParameterSet global, String key, ToolParameterSet... local ) throws IllegalValueException {
		Object o = global.getParameterForName(key).getValue();
		if( o != null ) {
			for( int i = 0; i < local.length; i++ ) {
				local[i].getParameterForName(key).setValue(o);
			}
		}
	}
	
	Thread killer;
	
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		stop=false;
		pipelineProtocol = 
				//new SysProtocol(); //XXX for debugging the Galaxy integration
				protocol;
		this.threads=threads;
		
		pipelineProtocol.append("Running the GeMoMaPipeline with " + threads + " threads\n\n" );
		
		//create temp dir
		File basic = new File(GeMoMa.GeMoMa_TEMP);
		File dir = Files.createTempDirectory(basic.toPath(), "GeMoMaPipeline-").toFile();
		dir.mkdirs();
		home = dir.toString() + "/";
		pipelineProtocol.append("temporary directory: " + home);
		
		Extractor extractor = new Extractor(maxSize);
		GeMoMa gemoma = new GeMoMa(maxSize,timeOut,maxTimeOut);
		GeMoMaAnnotationFilter gaf = new GeMoMaAnnotationFilter();

		ParameterSet pa = (ParameterSet) parameters.getParameterForName("RNA-seq evidence").getValue();
		rnaSeq=pa.getNumberOfParameters()>0;
		
		extractorParams = extractor.getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("Extractor parameter set")).getValue(), extractorParams);
		
		gemomaParams = gemoma.getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("GeMoMa parameter set")).getValue(), gemomaParams);
		
		gafParams = gaf.getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("GAF parameter set")).getValue(), gafParams);

		String[] key = {"selected", "genetic code"};
		for( int i = 0; i < key.length; i++ ) {
			setParameters(parameters, key[i], extractorParams, gemomaParams);
		}		
		setParameters(parameters, "tag", gemomaParams, gafParams );
		
		selected = Tools.getSelection( (String) parameters.getParameterForName(key[0]).getValue(), maxSize, protocol );
		
		eValue = (Double) gemomaParams.getParameterForName("e-value").getValue();
		gapOpen = (Integer) gemomaParams.getParameterForName("gap opening").getValue();
		gapExt = (Integer) gemomaParams.getParameterForName("gap extension").getValue();
		
		
		process = new ArrayList<Process>();
		list = new ArrayList<FlaggedRunnable>();
		queue = Executors.newFixedThreadPool(threads);
		ecs = new ExecutorCompletionService<>(queue);
		killer = new Thread() {
			public void run() {
				//pipelineProtocol.append("shut down hook\n");
				queue.shutdown();
                try {
                    if (!queue.awaitTermination(10, TimeUnit.MILLISECONDS)) {
                        //log.warn("Executor did not terminate in the specified time.");
                        //List<Runnable> droppedTasks = 
                    	queue.shutdownNow();
                        //log.warn("Executor was abruptly shut down. " + droppedTasks.size() + " tasks will not be executed.");
                    }
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }/**/
				for( int i = 0; i < process.size(); i++ ) {
					Process p = process.get(i);
					if( p.isAlive() ) {
						//System.out.println("shutdown thread: " + i);
						p.destroyForcibly();
					}
				}
			}
		};
		Runtime.getRuntime().addShutdownHook(killer);
		
		target = parameters.getParameterForName("target genome").getValue().toString();	
		
		//first part: preparation
		addNewPhase();
		
		//makeblastdb
		add(new JmakeBlastDB());
		
		//extract gene models from reference
		ExpandableParameterSet ref = (ExpandableParameterSet) ((ParameterSetContainer) parameters.getParameterForName("reference species")).getValue();
		int sp = ref.getNumberOfParameters();
		species = new ArrayList<Species>();
		for( int s=0; s < sp; s++){			
			SelectionParameter sel = (SelectionParameter) ((ParameterSetContainer)ref.getParameterAt(s)).getValue().getParameterAt(0);
			SimpleParameterSet currentRef = (SimpleParameterSet) sel.getValue();
			switch( currentRef.getNumberOfParameters() ) {
				case 1: //BUSCO
					String path = buscoPath+currentRef.getParameterAt(0).getValue() + "/";
					File clade = new File(path);
					String[] speciesDirs = clade.list(dirFilter);
					for( int i = 0; i < speciesDirs.length; i++ ) {
						add(new JExtractorAndSplit(
								speciesDirs[i],
								path + speciesDirs[i] + "/" + Extractor.name[0] + "." + Extractor.type[0], 
								path + speciesDirs[i] + "/" + Extractor.name[1] + "." + Extractor.type[1],
								null
						));
					}
					break;
				case 3: //from reference genome & annotation
					add(new JExtractorAndSplit(
							(String) currentRef.getParameterForName("ID").getValue(),
							currentRef.getParameterForName("annotation").getValue().toString(), 
							currentRef.getParameterForName("genome").getValue().toString()
					));
					break;
				case 4: //pre-extracted
					add(new JExtractorAndSplit(
							(String) currentRef.getParameterForName("ID").getValue(),
							currentRef.getParameterForName("cds parts").getValue().toString(), 
							(String) currentRef.getParameterForName("assignment").getValue(),
							(String) currentRef.getParameterForName("query proteins").getValue()
					));
					break;
				default: throw new UnsupportedOperationException();
			}
		}
		pipelineProtocol.append("Running GeMoMa for "+ species.size() + " reference species.\n");
		
		boolean nothing = false;
		res = new ArrayList<Result>();
		if( species.size() > 0 ) {
			//ERE
			rnaSeqData = new RNASeq();
			JEREAndFill ere;
			if( rnaSeq ) {
				Parameter p = pa.getParameterForName("introns");
				if( p != null ) {
					ExpandableParameterSet exp = (ExpandableParameterSet) p.getValue();
					for( int i = 0; i < exp.getNumberOfParameters(); i++ ) {
						rnaSeqData.introns.add( (String) ((SimpleParameterSet) exp.getParameterAt(i).getValue()).getParameterAt(0).getValue() );
					}
					
					exp = (ExpandableParameterSet) pa.getParameterForName("coverage").getValue();
					for( int i = 0; i < exp.getNumberOfParameters(); i++ ) {
						//SelectionParameter sel 
						SimpleParameterSet ps = (SimpleParameterSet) exp.getParameterAt(i).getValue();
						ps = (SimpleParameterSet) ps.getParameterAt(0).getValue();
						if( ps.getNumberOfParameters() == 1 ) {
							rnaSeqData.coverageUn.add( (String) ps.getParameterAt(0).getValue() );
						} else if( ps.getNumberOfParameters() == 2 ) {
							rnaSeqData.coverageFwd.add( (String) ps.getParameterAt(0).getValue() );
							rnaSeqData.coverageRC.add( (String) ps.getParameterAt(1).getValue() );
						}					
					}
					
					ere = new JEREAndFill();
				} else {
					ere = new JEREAndFill((ToolParameterSet) pa);
				}
			} else {
				ere = new JEREAndFill();
			}
			add(ere);
			
			//wait until first phase has been finished
			wait(list.size());
	
			
			//second phase = tblastn + GeMoMa
			if( !queue.isShutdown() ) {
				addNewPhase();
			
				//tblastn
				splitsPerSpecies = new IntList();
				RegExFilenameFilter filter = new RegExFilenameFilter("splits", Directory.FORBIDDEN, true, "split-.*fasta" );
				int anz=0;
				for( int s = 0; s < speciesCounter; s++ ) {
					File d = new File( home + s );
					String[] f = d.list(filter);
					splitsPerSpecies.add(f.length);
					pipelineProtocol.append("species " + s + ": " + f.length + " splits\n");
					for( int p = 0; p < f.length; p++ ) {
						add( new JTblastn(s, p) );
						anz++;
					}
				}
				
				//wait until second part has been finished
				wait(2*anz);
			}
			
			//third part = cat
			if( !queue.isShutdown() ) {
				addNewPhase();
				for( int s = 0; s < speciesCounter; s++ ) {
					add( new JCat(s) );
				}
				
				//wait until third part has been finished
				wait(species.size());/**/
			}
			
			//final prediction
			if( !queue.isShutdown() ) {
				addNewPhase();
				add(new JGAF());
				wait(1);
			}
			
			//TODO AnnotationFinalizer
			
			//extract predictions
			boolean cds = (Boolean) parameters.getParameterForName("predicted CDSs").getValue();
			boolean protein = (Boolean) parameters.getParameterForName("predicted proteins").getValue();
			boolean genomic = (Boolean) parameters.getParameterForName("predicted genomic regions").getValue();
			if( !queue.isShutdown() & (cds | protein | genomic) ) {
				addNewPhase();
				add(new JExtractor(cds, protein, genomic));
				wait(1);
			}
			
			//needs to be done
			if( !queue.isShutdown() ) {
				queue.shutdown();
				queue.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			}
		} else {
			wait(1);
			nothing=true;
		}
		
		//statistics
		HashMap<String,int[]> stat = new HashMap<String, int[]>();
		for( int i = 0; i < list.size(); i++ ) {
			FlaggedRunnable f = list.get(i);
			int[] val = stat.get(f.getClass().getName());
			if( val  == null ) {
				val = new int[JobStatus.values().length];
				stat.put(f.getClass().getName(), val);
			}
			val[f.status.ordinal()]++;
		}
		
		pipelineProtocol.append("\nStatistics:\n");
		pipelineProtocol.append("Job");
		JobStatus[] st = JobStatus.values();
		for( int i = 0; i < st.length; i++ ) {
			pipelineProtocol.append("\t" + st[i] );
		}
		pipelineProtocol.append("\n");
		Iterator<Entry<String,int[]>> it = stat.entrySet().iterator();
		int success = 0;
		while( it.hasNext() ) {
			Entry<String,int[]> e = it.next();
			int[] val = e.getValue();
			String name = e.getKey();
			int idx = name.lastIndexOf("$");
			pipelineProtocol.append( (idx>0 ? name.substring(idx+2) : name) );
			for( int i = 0; i < val.length; i++ ) {
				pipelineProtocol.append( "\t" + val[i] );
			}
			pipelineProtocol.append( "\n");
			success+=val[4];
		}
		pipelineProtocol.append("\n");
		
		ToolResult result;
		if( !nothing && success==list.size() ) {
			pipelineProtocol.append("No errors detected.\n");
			for( int i = 0; i < speciesCounter; i++ ) {
				String unfiltered = home+"/" + i + "/unfiltered-predictions.gff";
				//pipelineProtocol.append(unfiltered + "\t" + (new File(unfiltered)).exists() +"\n" );
				FileRepresentation fr = new FileRepresentation(unfiltered);
				fr.getContent();//this is important otherwise the output is null, as the files will be deleted within the next lines
				res.add( new TextResult("unfiltered predictions from species "+i, "Result", fr, "gff", gemoma.getToolName(), null, true) );
			}
			
			result = new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
		} else {
			pipelineProtocol.append( (list.size()-success) + " jobs did not finish as expected. Please check the output carefully.\n");
			result = null;
		}		
		
		
		if( result != null || (boolean) parameters.getParameterForName("debug").getValue() ) {
			//delete
			// for avoiding erroneously deleting important files: 2 tricks have been implemented:
			// a) temporary folder GeMoMa_TEMP
			// b) RegExFilenameFilter
			try {
				Files.walkFileTree(new File( home ).toPath(), new SimpleFileVisitor<Path>() {
					RegExFilenameFilter filter = new RegExFilenameFilter("only relevant files", Directory.ALLOWED, true, ".*fasta", ".*gff", ".*txt", "blastdb.*", ".*tabular", ".*bedgraph");
					
					@Override
					public FileVisitResult visitFile(Path file, BasicFileAttributes attrs) throws IOException {
						return delete( file );
					}
		
					@Override
					public FileVisitResult postVisitDirectory(Path dir, IOException exc) throws IOException {
						return delete(dir);
					}
					
					FileVisitResult delete( Path p ) throws IOException {
						if( filter.accept( p.toFile() ) )  {
							Files.delete(p);
						}
						return FileVisitResult.CONTINUE;
					}
				});
			} catch( IOException io ) {
				protocol.appendWarning("Could not delete all temporary files.");
				protocol.appendThrowable(io);
			}
		}
		
		if( result == null ) {
			pipelineProtocol.flush();
			//Thread.sleep(1000);
			throw new RuntimeException("Did not finish as intended");
		} else {
			return result;
		}
	}
	
	void add( FlaggedRunnable f ) {
		if( !queue.isShutdown() ) {
			list.add(f);
			ecs.submit(f,null);
		}
	}
	
	void wait( int num ) throws InterruptedException, ExecutionException {
		for (int i = 0; i < num; i++) {
			  ecs.take();
			  if( stop ) break;
		}
	}
	
	static class QuiteSysProtocol extends SysProtocol {
		public QuiteSysProtocol() {
			out = new PrintStream( SafeOutputStream.getSafeOutputStream(null) );
		}
	}
	
	/**
	 * Enum for the status of {@link FlaggedRunnable} jobs.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see FlaggedRunnable
	 */
	static enum JobStatus {
		WAITING,
		RUNNING,
		INTERRUPTED,
		FAILED,
		SUCCEEDED;
	}
	
	/**
	 * Abstract class for all jobs.
	 * 
	 * @author Jens Keilwagen
	 */
	abstract class FlaggedRunnable implements Runnable {
		JobStatus status;
		
		public FlaggedRunnable() {
			status = JobStatus.WAITING;
		}

		public final void run() {
			status = JobStatus.RUNNING;
			try {
				doJob();
			} catch (InterruptedException e ) {
				status = JobStatus.INTERRUPTED;
			} catch ( Throwable e ) {
				pipelineProtocol.appendThrowable( e );
				status = JobStatus.FAILED;
			}
			if( status == JobStatus.FAILED && !queue.isShutdown() ) {
				queue.shutdown();
				stop=true;
				killer.run();
			}
			if( status == JobStatus.RUNNING ) {
				status = JobStatus.SUCCEEDED;
			}
		}
		
		public abstract void doJob() throws Exception;
		
		int runProcess( String... cmd ) throws Exception {
			ProcessBuilder pb = new ProcessBuilder(cmd);
			pb.redirectErrorStream(true);//TODO ?
			Process p = pb.start();
			InputStream in = p.getInputStream();
			process.add(p);
			p.waitFor();
						
			if( in.available() > 0 ) {
				byte[] b = new byte[in.available()];
				in.read(b);
				pipelineProtocol.appendWarning(new String(b));
			}
			return p.exitValue();
		}
	}
	
	static String escape( String s ) {
		if( s.indexOf(' ') < 0 ) {
			return s;
		} else {
			return "\\\""+s+"\\\"";
		}
	}
	
	/**
	 * Job for running the ERE.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see ExtractRNAseqEvidence
	 */
	class JEREAndFill extends FlaggedRunnable {
		ToolParameterSet params;
		
		JEREAndFill() {
			params = null;			
		}
		
		JEREAndFill( ToolParameterSet ps ) {
			params = ps;			
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting ERE\n");
			Protocol protocol;
			if( params != null ) {
				protocol = new QuiteSysProtocol();
				//run ERE
				ExtractRNAseqEvidence ere = new ExtractRNAseqEvidence();
				ToolResult res = ere.run(params, protocol, new ProgressUpdater(), 1);
				CLI.writeToolResults(res, (SysProtocol) protocol, home+"/", ere, params);
				
				rnaSeqData.introns.add(home + "/introns.gff");
				
				boolean cov = ((Boolean) params.getParameterForName("coverage").getValue());
				if( cov ) {
					String fName = home + "/coverage.bedgraph";
					if( new File(fName).exists() ) {//unstranded
						rnaSeqData.coverageUn.add(fName);
					} else {//stranded
						rnaSeqData.coverageFwd.add(home + "/coverage_forward.bedgraph");
						rnaSeqData.coverageRC.add(home + "/coverage_reverse.bedgraph");
					}
				}
			}
			//prepare GeMoMa
			protocol = GeMoMaPipeline.this.pipelineProtocol;
			
			setGeMoMaParams(protocol);
			GeMoMa.fill( protocol, false, maxSize, 
					(String) gemomaParams.getParameterForName("target genome").getValue(), 
					(String) gemomaParams.getParameterForName("selected").getValue(),
					(Integer) gemomaParams.getParameterForName("reads").getValue(),
					(ExpandableParameterSet)((ParameterSetContainer)gemomaParams.getParameterAt(5)).getValue(),
					(ExpandableParameterSet)((ParameterSetContainer)gemomaParams.getParameterAt(8)).getValue()
			);			
			//protocol.append( "Fill okay\n" );
		}	
	}	

	/**
	 * Job for running makeblastdb.
	 * 
	 * @author Jens Keilwagen
	 */
	class JmakeBlastDB extends FlaggedRunnable {
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting makeblastdb\n");
			//log file is necessary as otherwise Java sometimes does not finish waitFor() 
			int exitCode = runProcess( "makeblastdb", "-out", escape(home+"/blastdb"), "-hash_index", "-in", escape(target), "-title", "target", "-dbtype", "nucl", "-logfile", home+"/blastdb-logfile" );
			
			//read logfile
			BufferedReader r = new BufferedReader( new FileReader(home+"/blastdb-logfile") );
			String line;
			StringBuffer log = new StringBuffer();
			HashMap<String,int[]> problems = new HashMap<String, int[]>();
			while( (line=r.readLine()) != null ) {
				if( line.startsWith("Error:")) {
					line = line.replaceAll("about ..%", "about xx%");
					
					int[] num = problems.get(line);
					if( num == null ) {
						num = new int[1];
						problems.put(line, num);
					}
					num[0]++;
				} else {
					log.append(line + "\n");
				}
			}
			r.close();
			
			//summarize errors
			Iterator<Entry<String, int[]>> it = problems.entrySet().iterator();
			if( it.hasNext() ) {
				log.append("\nproblems:\n");
				while( it.hasNext() ) {
					Entry<String,int[]> e = it.next();
					int num = e.getValue()[0];
					if( num == 1 ) {
						log.append( "once: " + e.getKey() + "\n");
					} else {
						log.append( num + " times: " + e.getKey() + "\n");
					}
				}
			}
			pipelineProtocol.append(log.toString());
			if( exitCode> 0 ) {
				status=JobStatus.FAILED;
			}
		}	
	}
	
	/**
	 * Job for running the Extractor and the FastaSplitter on one species.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see Extractor
	 * @see FastaSplitter
	 */
	class JExtractorAndSplit extends FlaggedRunnable {
		int speciesIndex;
		ToolParameterSet params;
		
		private JExtractorAndSplit( String specName ) {
			species.add( new Species(specName) );
			speciesIndex = speciesCounter++;
			params=null;
		}
		
		JExtractorAndSplit( String specName, String annotation, String genome ) throws IllegalValueException, CloneNotSupportedException {
			this(specName);
			params = extractorParams.clone();
			params.getParameterForName("annotation").setValue(annotation);
			params.getParameterForName("genome").setValue(genome);			
		}
		
		JExtractorAndSplit( String specName, String cds_parts, String assignment, String protein ) {
			this(specName);
			Species sp = species.get(speciesIndex);
			sp.cds=cds_parts;
			sp.assignment=assignment;
			sp.protein = protein;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting extractor for species " + speciesIndex+"\n");
			Species sp = species.get(speciesIndex);
			String outDir = home + "/" + speciesIndex + "/";
			File home = new File( outDir );
			home.mkdirs();
			if( params != null ) {
				SysProtocol protocol = new QuiteSysProtocol();
				Extractor extractor = new Extractor(maxSize);
				ToolResult res = extractor.run(params, protocol, new ProgressUpdater(), 1);
				CLI.writeToolResults(res, protocol, outDir, extractor, params);
				
				sp.cds = outDir + Extractor.name[0] + "." + Extractor.type[0];
				sp.assignment = outDir + Extractor.name[1] + "." + Extractor.type[1];
				String fName = outDir+Extractor.name[2] + "." + Extractor.type[2];
				if( new File(fName).exists() ) {
					sp.protein = fName;
				}
			} else if( selected != null ) {
				//find genes
				BufferedReader r = new BufferedReader( new FileReader( sp.assignment ) );
				String line;
				HashSet<String> genes = new HashSet<String>();
				while( (line=r.readLine())!= null ) {
					if( line.charAt(0) == '#' ) continue;
					String[] split = line.split("\t");
					if( selected.containsKey(split[1]) ) {
						genes.add(split[0]);
					}
				}
				r.close();
				
				//filter cds parts
				r = new BufferedReader( new FileReader( sp.cds ) );
				sp.cds = outDir + Extractor.name[0] + "." + Extractor.type[0];
				BufferedWriter w = new BufferedWriter( new FileWriter( sp.cds ) );
				boolean keep = false;
				while( (line=r.readLine())!= null ) {
					if( line.length()== 0 ) {
						continue;
					} else if( line.charAt(0) == '>' ) {
						int idx = line.lastIndexOf("_");
						if( idx < 0 ) {
							idx = line.length();
						}
						String id = line.substring(1,idx);
						keep = genes.contains(id);
					}

					if( keep ) {
						w.append(line);
						w.newLine();
					}
				}
				r.close();
				w.close();
			}
			FastaSplitter.main(new String[]{sp.cds,""+threads,"_",outDir});
		}	
	}
	
	/**
	 * Job for running a tblastn.
	 * 
	 * @author Jens Keilwagen
	 */
	class JTblastn extends FlaggedRunnable {
		int species, split;
		
		JTblastn( int species, int split ) throws IllegalValueException, CloneNotSupportedException {
			this.species=species;
			this.split = split;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting tblastn split=" + split + " for species " + species + "\n");
			int exitCode = runProcess( "tblastn", "-query", home+species+"/split-"+split+".fasta",
					"-db", escape(home+"blastdb"), "-evalue",""+eValue, "-out", home+species+"/tblastn-"+split+".txt",
					"-outfmt", "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles",
					"-db_gencode","1", "-matrix", "BLOSUM62", "-seg", "no", "-word_size", "3", "-comp_based_stats", "F", "-gapopen", ""+gapOpen, "-gapextend", ""+gapExt, "-num_threads", "1" );
			
			if( exitCode <= 1 ) { //Warnings/error in query will be ignored: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.Tc/
				if( !queue.isShutdown() ) {
					add(new JGeMoMa(species, split));
				}
			} else {
				status = JobStatus.FAILED;
			}
		}
	}
	
	/**
	 * Job for running GeMoMa on one split of one species.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see GeMoMa
	 */
	class JGeMoMa extends FlaggedRunnable {
		int species, split;
		
		JGeMoMa( int species, int split ) {
			this.species = species;
			this.split = split;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting GeMoMa split=" + split + " for species " + species + "\n");
			
			Species sp = GeMoMaPipeline.this.species.get(species);
			ToolParameterSet params = gemomaParams.clone();
			params.getParameterForName("tblastn results").setValue(home + "/" + species + "/tblastn-"+split+".txt");
			params.getParameterForName("cds parts").setValue(sp.cds);
			params.getParameterForName("assignment").setValue(sp.assignment);
			if( sp.protein != null ) {
				params.getParameterForName("query proteins").setValue(sp.protein);
			}
			
			SysProtocol protocol = new QuiteSysProtocol();
			GeMoMa gemoma = new GeMoMa(maxSize,timeOut,maxTimeOut);
			ToolResult res = gemoma.run(params, protocol, new ProgressUpdater(), 1);
			String outDir = home + "/" + species + "/"+split +"/";
			CLI.writeToolResults(res, protocol, outDir, gemoma, params);
		}
	}

	/**
	 * Job for concatenating predictions of one species
	 * 
	 * @author Jens Keilwagen
	 */
	class JCat extends FlaggedRunnable {
		int speciesIndex;
		
		JCat( int s ) {
			speciesIndex = s;
		}
		
		@Override
		public void doJob () throws Exception {
			pipelineProtocol.append("starting cat for species " + speciesIndex + "\n");
			BufferedWriter w= new BufferedWriter( new FileWriter(home+"/" + speciesIndex + "/unfiltered-predictions.gff") );
			for( int sp = 0; sp <splitsPerSpecies.get(speciesIndex); sp++ ) {
				BufferedReader r = new BufferedReader( new FileReader(home+"/" + speciesIndex + "/" + sp + "/predicted_annotation.gff") );
				String line;
				while( (line=r.readLine()) != null ) {
					if( line.charAt(0) != '#' || sp==0 ) {
						w.append(line);
						w.newLine();
					}
				}
				r.close();
			}
			w.close();
		}	
	}
	
	/**
	 * Job for finally combining the predictions
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see GAF
	 */
	class JGAF extends FlaggedRunnable {
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting GAF\n");
			SysProtocol protocol = new QuiteSysProtocol();
			
			GeMoMaAnnotationFilter gaf = new GeMoMaAnnotationFilter();
			
			ExpandableParameterSet eps = (ExpandableParameterSet) gafParams.getParameterForName("predicted annotation").getValue();
			int i = 0;
			for( int s = 0; s < speciesCounter; s++ ) {
				if( s>0 ) {
					eps.addParameterToSet();
				}
				ParameterSet ps = ((ParameterSetContainer) eps.getParameterAt(i++)).getValue();
				String id = species.get(s).id;
				if( id != null ) {
					ps.getParameterForName("prefix").setValue(id);
				}
				ps.getParameterForName("gene annotation file").setValue(home+"/" + s + "/unfiltered-predictions.gff");
			}
			
			ToolResult tr = gaf.run(gafParams, protocol, new ProgressUpdater(), 1);
			res.add( tr.getRawResult()[0].getResultAt(0) );
			CLI.writeToolResults(tr, protocol, home, gaf, gafParams);
		}	
	}
	
	/**
	 * Job for running the Extractor on the final prediction.
	 * 
	 * @author Jens Keilwagen
	 * 
	 * @see Extractor
	 */
	class JExtractor extends FlaggedRunnable {
		Extractor extractor;
		ToolParameterSet params;
				
		JExtractor( boolean cds, boolean protein, boolean genomic ) throws IllegalValueException, CloneNotSupportedException {
			extractor = new Extractor(-1);
			params = extractor.getToolParameters();
			
			params.getParameterForName("annotation").setValue(home + "/filtered_predictions.gff");
			params.getParameterForName("genome").setValue(target);
			params.getParameterForName("Ambiguity").setValue("AMBIGUOUS");
			
			params.getParameterForName(Extractor.name[2]).setValue(protein);
			params.getParameterForName(Extractor.name[3]).setValue(cds);
			params.getParameterForName(Extractor.name[4]).setValue(genomic);
		}
				
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting extractor for final prediction\n");
			SysProtocol protocol = new QuiteSysProtocol();
			ToolResult tr = extractor.run(params, protocol, new ProgressUpdater(), 1);
			ResultSet raw = tr.getRawResult()[0];
			for( int i = 2; i < raw.getNumberOfResults(); i++ ) {
				res.add( raw.getResultAt(i) );
			}
			CLI.writeToolResults(tr, protocol, home, extractor, params);
		}	
	}

	@Override
	public String getToolName() {
		return "GeMoMa pipeline";
	}

	@Override
	public String getToolVersion() {
		return GeMoMa.VERSION;
	}

	@Override
	public String getShortName() {
		return "GeMoMaPipeline";
	}

	@Override
	public String getDescription() {
		return "runs the complete GeMoMa pipeline for one target species and multiple refernece species";
	}

	@Override
	public String getHelpText() {
		return "**What it does**\n\nThis tool can be used to run the complete GeMoMa pipeline. The tool is multi-threaded and can utilize all compute cores on one machine, but not distributed as for instance in a compute cluster.\nIt basically runs: makeblastdb, Extract RNA-seq evidence (ERE), Extractor, FastaSplitter, tblastn, GeMoMa, cat, and GeMoMa Annotation Filter (GAF).\n\n"
				+ GeMoMa.REF;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", "filtered predictions"),
		};
	}
}