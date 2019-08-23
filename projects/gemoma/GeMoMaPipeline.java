/*
 * This file is part of Jstacs.
 * 
 * Jstacs is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Jstacs is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Jstacs. If not, see <http://www.gnu.org/licenses/>.
 * 
 * For more information on Jstacs, visit http://www.jstacs.de
 */

package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.SimpleFileVisitor;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
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
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.tools.ui.cli.CLI;
import de.jstacs.tools.ui.cli.CLI.QuietSysProtocol;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;
import de.jstacs.tools.ui.galaxy.Galaxy;
import de.jstacs.utils.Time;
import projects.FastaSplitter;

/**
 * The GeMoMa pipeline as one tool.
 * 
 * @author Jens Keilwagen
 * 
 * @see ExtractRNAseqEvidence
 * @see Denoise
 * @see Extractor
 * @see GeMoMa
 * @see GAF
 * @see AnnotationFinalizer
 */
public class GeMoMaPipeline extends GeMoMaModule {

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
			maxSize = -1;
			timeOut = 3600;
			maxTimeOut = 60*60*24*7;
		}
	}

	String home;
	String target;
	int threads, gapOpen, gapExt;
	double eValue;
	boolean rnaSeq, stop, all;
	
	ToolParameterSet params, denoiseParams, extractorParams, gemomaParams, gafParams, afParams;
	ExecutorService queue;
	ExecutorCompletionService ecs;
	
	ArrayList<Species> species;
	int speciesCounter = 0;
	
	RNASeq rnaSeqData;
	
	ArrayList<Process> process;
	ArrayList<FlaggedRunnable> list;
	
	HashMap<String,String> selected;
	
	public ParameterSet getRelevantParameters( ParameterSet params, String... remove ) {
		HashSet<String> removeNames = new HashSet<String>();
		for( String r : remove ) {
			removeNames.add(r);
		}
		ArrayList<Parameter> list = new ArrayList<Parameter>();
		for( int i = 0; i < params.getNumberOfParameters(); i++ ) {
			Parameter p = params.getParameterAt(i);
			if( !removeNames.contains(p.getName()) ) {
				if( p instanceof SelectionParameter ) {
					SelectionParameter sel = (SelectionParameter) p;
					int selected = sel.getSelected();
					ParameterSet ps = sel.getParametersInCollection();
					try {
						for( int j = 0; j < ps.getNumberOfParameters(); j++ ) {
							Parameter q = ps.getParameterAt(j);
							if( q instanceof ParameterSetContainer ) {
								ParameterSetContainer psc = (ParameterSetContainer) q;
								psc.setValue( getRelevantParameters( psc.getValue(), remove ) );
								sel.setValue(q);
							}							
						}
						sel.setValue(ps.getParameterAt(selected).getName());
					} catch (IllegalValueException e) {
						System.out.println(p.getName());
						e.printStackTrace();
						System.exit(1);
					}

/*					ParameterSet old = (ParameterSet) p.getValue();
					ParameterSet ps = getRelevantParameters(old, remove);
					try {
						p.setValue( ps );
					} catch (IllegalValueException e) {
						System.out.println(p.getName());
						
						System.out.println(ParameterSet.getName(old));
						System.out.println(ParameterSet.getName(ps));
						
						System.out.println(old.getNumberOfParameters());
						System.out.println(ps.getNumberOfParameters());
						e.printStackTrace();
						System.exit(1);
					}*/
				}
				list.add(p);
			}			
		}
		if( params instanceof ToolParameterSet ) {
			ToolParameterSet tps = (ToolParameterSet) params;
			return new ToolParameterSet( tps.getToolName(), list.toArray(new Parameter[0])  );
		} else {
			return new SimpleParameterSet( list.toArray(new Parameter[0]) );
		}
	}
	
	private static String buscoPath= "BUSCO-references/";
	private static FilenameFilter dirFilter = new RegExFilenameFilter("", Directory.REQUIRED, true, ".*" );
	
	public ToolParameterSet getToolParameters() {
		ToolParameterSet ere = new ExtractRNAseqEvidence().getToolParameters();
		ParameterSet denoise = getRelevantParameters(new Denoise().getToolParameters(), "introns", "coverage", "maximum intron length" );
		ParameterSet ex = getRelevantParameters(new Extractor(maxSize).getToolParameters(), "annotation", "genome", "selected", "verbose", "genetic code", Extractor.name[3], Extractor.name[4] );
		ParameterSet gem = getRelevantParameters(new GeMoMa(maxSize,timeOut,maxTimeOut).getToolParameters(), "search results", "target genome", "cds parts", "assignment", "query proteins", "selected", "verbose", "genetic code", "tag", "coverage", "introns", "sort", "maximum intron length" );
		ParameterSet gaf = getRelevantParameters(new GeMoMaAnnotationFilter().getToolParameters(), "predicted annotation", "tag");
		ParameterSet af = getRelevantParameters(new AnnotationFinalizer().getToolParameters(), "genome", "annotation", "tag", "introns", "reads", "coverage" );
		//XXX ex.getParameterForName("proteins").setDefault(true);
		ArrayList<String> keys = new ArrayList<String>();
		keys.add("own");
		keys.add("pre-extracted");
		
		try{
			ArrayList<SimpleParameterSet> values = new ArrayList<SimpleParameterSet>();
			values.add( new SimpleParameterSet(
							new SimpleParameter(DataType.STRING,"ID","ID to distinguish the different reference species", false, ""),
							new FileParameter( "annotation", "Reference annotation file (GFF or GTF), which contains gene models annotated in the reference genome", "gff,gtf", true ),
							new FileParameter( "genome", "Reference genome file (FASTA)", "fasta,fa.gz,fasta.gz,fa.gz",  true ),
							new SimpleParameter(DataType.DOUBLE,"weight","the weight can be used to prioritize predictions from different input files; each prediction will get an additional attribute sumWeight that can be used in the filter", false, new NumberValidator<Double>(0d, 1000d), 1d)
			) );
			values.add( new SimpleParameterSet(
							new SimpleParameter(DataType.STRING,"ID","ID to distinguish the different reference species", false, ""),
							new FileParameter( "cds parts", "The query cds parts file (FASTA), i.e., the cds parts that have been blasted", "fasta", true ),
							new FileParameter( "assignment", "The assignment file, which combines parts of the CDS to transcripts", "tabular", false ),
							new FileParameter( "query proteins", "optional query protein file (FASTA) for computing the optimal alignment score against complete protein prediction", "fasta", false ),
							new SimpleParameter(DataType.DOUBLE,"weight","the weight can be used to prioritize predictions from different input files; each prediction will get an additional attribute sumWeight that can be used in the filter", false, new NumberValidator<Double>(0d, 1000d), 1d)
			) );
			
			/*TODO BUSCO
			File busco = new File( buscoPath );
			if( busco.exists() && busco.isDirectory() ) {
				File[] d = busco.listFiles(dirFilter);
				
				String[] k = new String[d.length];
				for( int i = 0; i < d.length; i++ ) {
					k[i] = d[i].getName();
				}
				keys.add("busco");
				values.add( new SimpleParameterSet( new SelectionParameter(DataType.STRING, k, k, "BUSCO clades", "the clade to be used", true) ) );
			}/**/
			
			return new ToolParameterSet( getShortName(),
				new FileParameter( "target genome", "Target genome file (FASTA)", "fasta,fa.gz,fasta.gz,fa.gz",  true ),

				new ParameterSetContainer( "reference species", "", 
						new ExpandableParameterSet( new SimpleParameterSet(	
								new SelectionParameter( DataType.PARAMETERSET, keys.toArray(new String[0]),
										values.toArray(new SimpleParameterSet[0]), "species", "data for reference species", true
								)
						), "reference", "", 1 )
				),/**/
				
				new FileParameter( "selected", "The path to list file, which allows to make only a predictions for the contained transcript ids. The first column should contain transcript IDs as given in the annotation. Remaining columns can be used to determine a target region that should be overlapped by the prediction, if columns 2 to 5 contain chromosome, strand, start and end of region", "tabular,txt", maxSize>-1 ),
				new FileParameter( "genetic code", "optional user-specified genetic code", "tabular", false ),
				new SimpleParameter( DataType.INT, "maximum intron length", "The maximum length of an intron", true, 15000 ),
				new SimpleParameter( DataType.BOOLEAN, "tblastn", "if *true* tblastn is used as search algorithm, otherwise mmseqs is used. Tblastn and mmseqs need to be installed to use the corresponding option", true, true),
				new SimpleParameter( DataType.STRING, "tag", "A user-specified tag for transcript predictions in the third column of the returned gff. It might be beneficial to set this to a specific value for some genome browsers.", true, GeMoMa.TAG ),
				
				new SelectionParameter( DataType.PARAMETERSET, 
						new String[]{"NO", "MAPPED","EXTRACTED"},
						new Object[]{
								new SimpleParameterSet(),
								ere,
								new SimpleParameterSet(
										new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
												new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff", true )
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
				new SelectionParameter( DataType.PARAMETERSET, 
						new String[]{"DENOISE", "RAW"},
						new Object[]{
								denoise,
								new SimpleParameterSet()
						},
						"denoise", "removing questionable introns that have been extracted by ERE", true ),
				
				new ParameterSetContainer("Extractor parameter set", "parameters for the Extrator module of GeMoMa", ex ),
				new ParameterSetContainer("GeMoMa parameter set", "parameters for the GeMoMa", gem ),
				new ParameterSetContainer("GAF parameter set", "parameters for the GAF module of GeMoMa", gaf ),
				new ParameterSetContainer("AnnotationFinalizer parameter set", "parameters for the AnnotationFinalizer module of GeMoMa", af ),
				
				new SimpleParameter( DataType.BOOLEAN, "predicted proteins", "If *true*, returns the predicted proteins of the target organism as fastA file", true, true ),
				new SimpleParameter( DataType.BOOLEAN, "predicted CDSs", "If *true*, returns the predicted CDSs of the target organism as fastA file", true, false ),
				new SimpleParameter( DataType.BOOLEAN, "predicted genomic regions", "If *true*, returns the genomic regions of predicted gene models of the target organism as fastA file", true, false ),
				
				new SimpleParameter( DataType.BOOLEAN, "output individual predictions", "If *true*, returns the predictions for each reference species", true, false ),
				
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
				
				if( b instanceof SelectionParameter ) {
					int idx = ((SelectionParameter)a).getSelected();
					b.setValue( a.getParametersInCollection().getParameterAt(idx) );
				}
				
				setParameters( ps, (ParameterSet) b.getValue() );
			} else {
				q.setValue( p.getValue() );
			}
		}
	}
	
	void setGeMoMaParams( Protocol protocol ) throws IllegalValueException, CloneNotSupportedException {
		gemomaParams.getParameterForName("target genome").setValue(target);
		ParameterSet ps = (ParameterSet) params.getParameterForName("GeMoMa parameter set").getValue();
		gemomaParams.getParameterForName("reads").setValue( ps.getParameterForName("reads").getValue() );
		setRNASeqParams( gemomaParams, protocol );
	}
	
	void setRNASeqParams( ToolParameterSet p, Protocol protocol ) throws IllegalValueException, CloneNotSupportedException {
		//introns
		ExpandableParameterSet exp = (ExpandableParameterSet) ((ParameterSetContainer) p.getParameterForName("introns")).getValue();
		for( int i = 0; i < rnaSeqData.introns.size(); i++ ) {
			if( exp.getNumberOfParameters() <= i ) {
				exp.addParameterToSet();
			}
			String current = rnaSeqData.introns.get(i);
			if( current != null ) {
				((ParameterSetContainer) exp.getParameterAt(i)).getValue().getParameterAt(0).setValue(current);
			}
		}

		//coverage
		exp = (ExpandableParameterSet) ((ParameterSetContainer) p.getParameterForName("coverage")).getValue();
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
		Double weight;
		String[] cds_parts;
		String[] searchResults;
		boolean hasCDS;
		
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
				if( local[i] != null ) {
					local[i].getParameterForName(key).setValue(o);
				}
			}
		}
	}
	
	Thread killer;

	static String mmseqs;
	
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		Thread.sleep(1000);
		Time t = Time.getTimeInstance(null);

		params=parameters;
		all=(Boolean) params.getParameterForName("output individual predictions").getValue();
		
		stop=false;
		pipelineProtocol = 
				//new SysProtocol(); //XXX for debugging the Galaxy integration
				protocol;
		this.threads=threads;
		
		String os = System.getProperty("os.name");
		if( os.startsWith("Windows") ) {
			mmseqs = "mmseqs.bat";
		} else {
			mmseqs = "mmseqs";
		}
		
		pipelineProtocol.append("Running the GeMoMaPipeline with " + threads + " threads\n\n" );
		
		//create temp dir
		File basic = new File(Tools.GeMoMa_TEMP);
		File dir = Files.createTempDirectory(basic.toPath(), "GeMoMaPipeline-").toFile();
		dir.mkdirs();
		home = dir.toString() + "/";
		pipelineProtocol.append("temporary directory: " + home + "\n\n");
		
		//checking whether the search algorithm software is installed
		boolean blast=(Boolean) parameters.getParameterForName("tblastn").getValue();
		String[] cmd;
		if( blast ) {
			cmd=new String[] {"tblastn","-version"};
		} else {
			cmd=new String[] {mmseqs,"version"};
		}
		ProcessBuilder pb = new ProcessBuilder(cmd);
		pb.redirectErrorStream(true);
		
		Process pr = pb.start();
        BufferedReader br = new BufferedReader(new InputStreamReader(pr.getInputStream()));
        String line = null;
        pipelineProtocol.append("search algorithm:\n");
        while ( (line = br.readLine()) != null) {
        	pipelineProtocol.appendWarning( "[" + cmd[0] + "]: " + line+"\n");
		}
        br.close();
		pr.waitFor();
		
		if( pr.exitValue() != 0 ) {
			System.out.println("Search algorithm not available!");
			System.exit(1);
		}

		Denoise denoise = new Denoise();		
		Extractor extractor = new Extractor(maxSize);
		GeMoMa gemoma = new GeMoMa(maxSize,timeOut,maxTimeOut);

		ParameterSet pa = (ParameterSet) parameters.getParameterForName("RNA-seq evidence").getValue();
		rnaSeq=pa.getNumberOfParameters()>0;

		ParameterSet dps = (ParameterSet) parameters.getParameterForName("denoise").getValue();
		if( rnaSeq && dps.getNumberOfParameters()>0 ) {
			denoiseParams = denoise.getToolParameters();
			setParameters(dps, denoiseParams);
		} else {
			denoiseParams = null;
		}
		
		extractorParams = extractor.getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("Extractor parameter set")).getValue(), extractorParams);
		
		gemomaParams = gemoma.getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("GeMoMa parameter set")).getValue(), gemomaParams);
		
		gafParams = new GeMoMaAnnotationFilter().getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("GAF parameter set")).getValue(), gafParams);

		afParams = new AnnotationFinalizer().getToolParameters();
		setParameters(((ParameterSetContainer) parameters.getParameterForName("AnnotationFinalizer parameter set")).getValue(), afParams);
		
		String[] key = {"selected", "genetic code"};
		for( int i = 0; i < key.length; i++ ) {
			setParameters(parameters, key[i], extractorParams, gemomaParams);
		}		
		setParameters(parameters, "tag", gemomaParams, gafParams, afParams );
		setParameters(parameters, "maximum intron length", gemomaParams, denoiseParams );
		
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
	
		//create target db
		if( blast ) {
			add(new JMakeBlastDB());
		} else {
			add( new JMmseqsCreateDB() );
		}
		
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
								null,
								path + speciesDirs[i] + "/" + Extractor.name[0] + "." + Extractor.type[0], 
								path + speciesDirs[i] + "/" + Extractor.name[1] + "." + Extractor.type[1],
								null,
								blast
						));
					}
					break;
				case 4: //from reference genome & annotation
					add(new JExtractorAndSplit(
							(String) currentRef.getParameterForName("ID").getValue(),
							(Double) currentRef.getParameterForName("weight").getValue(),
							currentRef.getParameterForName("annotation").getValue().toString(), 
							currentRef.getParameterForName("genome").getValue().toString(),
							blast
					));
					break;
				case 5: //pre-extracted
					add(new JExtractorAndSplit(
							(String) currentRef.getParameterForName("ID").getValue(),
							(Double) currentRef.getParameterForName("weight").getValue(),
							(String) currentRef.getParameterForName("cds parts").getValue().toString(), 
							(String) currentRef.getParameterForName("assignment").getValue(),
							(String) currentRef.getParameterForName("query proteins").getValue(),
							blast
					));
					break;
				default: throw new UnsupportedOperationException();
			}
		}
		//pipelineProtocol.append("Running GeMoMa for "+ species.size() + " reference species.\n");
		
		int usedSpec = 0;
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
			waitPhase();
	
			for( int s = 0; s < speciesCounter; s++ ) {
				Species current  = species.get(s);
				if( current.hasCDS ) {
					usedSpec++;
				}
			}
			if( usedSpec>0 ) {
				//second phase = tblastn + GeMoMa
				if( !queue.isShutdown() ) {
					addNewPhase();
					for( int s = 0; s < speciesCounter; s++ ) {
						Species current  = species.get(s);
						if( current.hasCDS ) {
							if( blast ) { //tblastn
								pipelineProtocol.append("species " + s + ": " + current.cds_parts.length + " splits\n");
								for( int p = 0; p < current.cds_parts.length; p++ ) {
									add( new JTblastn(s, p) );
								}
							} else { //mmseqs
								addMmseqsAndDummies(s);
							}
						}
					}
					//wait until second part has been finished
					waitPhase();
				}
				
				//third part = cat
				if( !queue.isShutdown() ) {
					addNewPhase();
					for( int s = 0; s < speciesCounter; s++ ) {
						Species current  = species.get(s);
						if( current.hasCDS ) {
							add( new JCat(s) );
						}
					}
					
					//wait until third part has been finished
					waitPhase();/**/
				}
				
				//filtered prediction
				if( !queue.isShutdown() ) {
					addNewPhase();
					add(new JGAF());
					waitPhase();
				}
				
				//UTR prediction
				if( !queue.isShutdown() ) {
					addNewPhase();
					add(new JAnnotationFinalizer());
					waitPhase();
				}
				
				//extract predictions
				boolean cds = (Boolean) parameters.getParameterForName("predicted CDSs").getValue();
				boolean protein = (Boolean) parameters.getParameterForName("predicted proteins").getValue();
				boolean genomic = (Boolean) parameters.getParameterForName("predicted genomic regions").getValue();
				if( !queue.isShutdown() & (cds | protein | genomic) ) {
					addNewPhase();
					add(new JExtractor(cds, protein, genomic));
					waitPhase();
				}
			}
			//needs to be done
			if( !queue.isShutdown() ) {
				queue.shutdown();
				queue.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
			}
		} else {
			wait(1);
		}
		
		//statistics
		HashMap<String,int[]> stat = new HashMap<String, int[]>();
		int anz = 0;
		for( int i = 0; i < list.size(); i++ ) {
			FlaggedRunnable f = list.get(i);
			if( !(f instanceof JSleep) ) {
				int[] val = stat.get(f.getClass().getName());
				if( val  == null ) {
					val = new int[JobStatus.values().length];
					stat.put(f.getClass().getName(), val);
				}
				val[f.status.ordinal()]++;
				anz++;
			}
		}
		
		pipelineProtocol.append("\nStatistics:\n");
		pipelineProtocol.append("Job");
		JobStatus[] st = JobStatus.values();
		for( int i = 0; i < st.length; i++ ) {
			pipelineProtocol.append("\t" + st[i] );
		}
		pipelineProtocol.append("\n");
		pipelineProtocol.append("---------------------------------------------------------\n");
		String[] keys = stat.keySet().toArray(new String[0]);
				
		Arrays.sort(keys, JobComparator.DEFAULT);
		int success = 0;
		for( int j = 0; j < keys.length; j++ ) {
			String name = keys[j];
			int[] val = stat.get(name);
			int idx = name.lastIndexOf("$");
			name = idx>0 ? name.substring(idx+2) : name;
			pipelineProtocol.append( /*JobComparator.getOrdinal(name) + "\t" +*/ name );
			for( int i = 0; i < val.length; i++ ) {
				pipelineProtocol.append( "\t" + val[i] );
			}
			pipelineProtocol.append( "\n");
			success+=val[4];
			
			if( j==1 || j==6 ) {
				pipelineProtocol.append( "\n");
			}
		}
		pipelineProtocol.append("\n");
		
		ToolResult result = null;
		if( usedSpec>0 && success==anz ) {
			pipelineProtocol.append("No errors detected.\n");
			if( all ) {
				for( int i = 0; i < speciesCounter; i++ ) {
					Species current = species.get(i);
					if( current.hasCDS ) {
						String unfiltered = home+"/" + i + "/unfiltered-predictions.gff";
						pipelineProtocol.append(unfiltered + "\t" + (new File(unfiltered)).exists() +"\n" );
						
						FileRepresentation fr = new FileRepresentation(unfiltered);
						//fr.getContent();//this is important otherwise the output is null, as the files will be deleted within the next lines
						//fr.setFilename(null);
						res.add( new TextResult("unfiltered predictions from species "+i, "Result", fr, "gff", gemoma.getToolName(), null, true) );
					}
				}
			}
			result = new ToolResult("", "", null, new ResultSet(res), parameters, getToolName(), new Date());
		} else {
			pipelineProtocol.appendWarning( (anz-success) + " jobs did not finish as expected. Please check the output carefully.\n");
		}		
		
		
		for( int i = 0; i < res.size(); i++ ) {
			Result r = res.get(i);
			if( r instanceof TextResult ) {
				TextResult tr = (TextResult) r;
				FileRepresentation fr = tr.getValue();
				fr.getContent();//this is important otherwise the output is null, as the files will be deleted within the next lines
				fr.setFilename(null);
			}
		}
		
		boolean debug = (boolean) parameters.getParameterForName("debug").getValue();
		if( result != null || !debug ) {
			//delete
			// for avoiding erroneously deleting important files: 2 tricks have been implemented:
			// a) temporary folder GeMoMa_TEMP
			// b) RegExFilenameFilter
			ArrayList<IOException> io = new ArrayList<IOException>();
			Files.walkFileTree(new File( home ).toPath(), new SimpleFileVisitor<Path>() {
				RegExFilenameFilter filter = new RegExFilenameFilter("only relevant files", Directory.ALLOWED, true, 
						".*fasta", ".*gff", ".*txt", ".*tabular", ".*bedgraph", ".*gff3", 
						"blastdb.*",
						"mmseqsdb.*", "t_orfs.*", "pref_.*", "aln.*", "translated_search\\.sh", "blastp\\.sh", "latest"
				);
				
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
						try {
							Files.delete(p);
						} catch( IOException e ) {
							io.add(e);
						}
					}
					return FileVisitResult.CONTINUE;
				}
			});
			if( io.size()>0 ) {
				protocol.appendWarning("Could not delete all temporary files.\n\n");
				for( int i = 0; i < io.size(); i++ ) {
					protocol.appendWarning(io.get(i).getMessage() + "\n");
				}
			}
		} else {
			protocol.appendWarning("Did not delete temporary files allowing to debug.\n\n");
		}
		
		long ti = Math.round(t.getElapsedTime());
		protocol.append("Elapsed time: " + ti + " seconds");
		int s = (int) (ti % 60);
		ti /= 60;
		int m = (int) (ti % 60);
		ti /= 60;
		protocol.append("\t(" + ti + "h " + m + "m " + s + "s)\n");
		
		if( result == null ) {
			pipelineProtocol.flush();
			//Thread.sleep(1000);
			RuntimeException r;
			if( tr == null ) {
				 r = new RuntimeException("Did not finish as intended." + (usedSpec==0?" No gene model was extracted from the references.":""));
			} else {
				r = new RuntimeException("Did not finish as intended. " + tr.getClass().getName() + ": " + tr.getMessage() );
				r.setStackTrace( tr.getStackTrace() );
			}
			throw r;
		} else {
			return result;
		}
	}
	
	static class JobComparator implements Comparator<String> {

		static JobComparator DEFAULT = new JobComparator();
		
		static HashMap<String,Integer> ordinal; 
		
		private JobComparator() {
			ordinal = new HashMap<String, Integer>();
			
			ordinal.put( "MakeBlastDB", 0);
			ordinal.put( "MmseqsCreateDB", 0 );
			ordinal.put( "EREAndFill", 1);
			
			ordinal.put( "ExtractorAndSplit", 2);
			ordinal.put( "Tblastn", 3);
			ordinal.put( "Mmseqs", 3);
			ordinal.put( "GeMoMa", 4);
			ordinal.put( "Cat", 5);
			ordinal.put( "GAF", 6);
			
			ordinal.put( "AnnotationFinalizer", 7);
			ordinal.put( "Extractor", 8);
		}
		
		static int getOrdinal( String s ) {
			int idx = s.lastIndexOf("$");
			s = idx>0 ? s.substring(idx+2) : s;
			Integer i = ordinal.get(s);
			return i==null? Integer.MAX_VALUE : i;
		}
		
		@Override
		public int compare(String o1, String o2) {
			return Integer.compare( getOrdinal(o1), getOrdinal(o2) );
		}
		
		
		
	}
	
	void add( FlaggedRunnable f ) {
		if( !queue.isShutdown() ) {
			list.add(f);
			ecs.submit(f,null);
		}
	}
	
	int finished = 0;
	
	void waitPhase() throws InterruptedException, ExecutionException {
		while( finished < list.size() ) {
//System.out.println("wait " + finished + " vs. " + list.size());
			  ecs.take();
			  if( stop ) break;
			  finished++;
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
	
	Throwable tr;
	
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
				tr=e;
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
			/*for( int i = 0; i < cmd.length; i++ ) {
				System.out.print( cmd[i] + " " );
			}
			System.out.println();*/
			
			ProcessBuilder pb = new ProcessBuilder(cmd);
			pb.redirectErrorStream(true);//TODO ?
			
			Process p = pb.start();
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String line = null;
            while ( (line = br.readLine()) != null) {
            	pipelineProtocol.appendWarning( "[" + cmd[0] + "]: " + line+"\n");
			}
            br.close();
			
			process.add(p);
			p.waitFor();
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
			Protocol protocol = new QuietSysProtocol();
			if( params != null ) {
				//run ERE
				ExtractRNAseqEvidence ere = new ExtractRNAseqEvidence();
				pipelineProtocol.append("starting ERE\n");
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
			
			if( denoiseParams != null ) {
				setRNASeqParams( denoiseParams, protocol );
				Denoise denoise = new Denoise();
				pipelineProtocol.append("starting Denoise\n");
				ToolResult res = denoise.run(denoiseParams, protocol, new ProgressUpdater(), 1);
				CLI.writeToolResults(res, (SysProtocol) protocol, home+"/", denoise, denoiseParams);
				
				rnaSeqData.introns.clear();
				rnaSeqData.introns.add(home + "/denoised_introns.gff");
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
	class JMakeBlastDB extends FlaggedRunnable {
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
	 * Job for running the mmseqs index on the target genome
	 * 
	 * @author Jens Keilwagen
	 */
	class JMmseqsCreateDB extends FlaggedRunnable {

		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting mmseqs createDB\n");
			int exitCode = runProcess( mmseqs, "createdb", escape(target), escape(home+"/mmseqsdb"), "--dont-split-seq-by-len", /*"--dont-shuffle", "--dbtype", "2",*/ "-v", "2" );
			
			if( exitCode> 0 ) {
				status=JobStatus.FAILED;
			}
			
			//pipelineProtocol.append( "Mmseqs createDB okay\n" );
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
		boolean split;
		
		private JExtractorAndSplit( String specName, Double w, boolean split ) {
			species.add( new Species(specName) );
			speciesIndex = speciesCounter++;
			species.get( speciesIndex ).weight = w;
			params=null;
			this.split=split;
		}
		
		JExtractorAndSplit( String specName, Double w, String annotation, String genome, boolean split ) throws IllegalValueException, CloneNotSupportedException {
			this(specName, w, split);
			params = extractorParams.clone();
			params.getParameterForName("annotation").setValue(annotation);
			params.getParameterForName("genome").setValue(genome);			
		}
		
		JExtractorAndSplit( String specName, Double w, String cds_parts, String assignment, String protein, boolean split ) {
			this(specName, w, split);
			Species sp = species.get(speciesIndex);
			sp.cds=cds_parts;
			sp.assignment=assignment;
			sp.protein=protein;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting extractor for species " + speciesIndex+"\n");
			Species sp = species.get(speciesIndex);
			String outDir = home + "/" + speciesIndex + "/";
			File home = new File( outDir );
			home.mkdirs();
			if( params != null ) {
				SysProtocol protocol = new QuietSysProtocol();
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
			File f = new File( sp.cds );
			sp.hasCDS = f.exists() && f.length()>0;
			if( sp.hasCDS ) {
				if( split && sp.hasCDS ) {
					FastaSplitter.main(new String[]{sp.cds,""+threads,"_",outDir});
					
					//Set files
					FilenameFilter filter = new RegExFilenameFilter("splits", Directory.FORBIDDEN, true, "split-.*fasta" );
					File d = new File( outDir );
					File[] h = d.listFiles(filter);
					sp.cds_parts = new String[h.length];
					for( int i = 0; i < h.length; i++ ) {
						sp.cds_parts[i] = h[i].getAbsolutePath();
					}
	
					sp.searchResults = new String[sp.cds_parts.length];
				}
			} else {
				pipelineProtocol.append("Did not extract any gene model from species " + speciesIndex + "\n");
			}
			//pipelineProtocol.append("Extract end\n");
		}	
	}
	
	/**
	 * Job for running tblastn.
	 * 
	 * @author Jens Keilwagen
	 */
	class JTblastn extends FlaggedRunnable {
		int speciesIndex, split;
		
		JTblastn( int species, int split ) throws IllegalValueException, CloneNotSupportedException {
			this.speciesIndex=species;
			this.split = split;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting tblastn split=" + split + " for species " + speciesIndex + "\n");
			Species current = species.get(speciesIndex);
			int exitCode = runProcess( "tblastn", "-query", current.cds_parts[split],
					"-db", escape(home+"blastdb"), "-evalue",""+eValue, "-out", home+speciesIndex+"/tblastn-"+split+".txt",
					"-outfmt", "6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles",
					"-db_gencode","1", "-matrix", "BLOSUM62", "-seg", "no", "-word_size", "3", "-comp_based_stats", "F", "-gapopen", ""+gapOpen, "-gapextend", ""+gapExt, "-num_threads", "1" );
			current.searchResults[split] = home+speciesIndex+"/tblastn-"+split+".txt";
						
			if( exitCode <= 1 ) { //Warnings/error in query will be ignored: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.Tc/
				if( !queue.isShutdown() ) {
					add(new JGeMoMa(speciesIndex, split));
				}
			} else {
				status = JobStatus.FAILED;
			}
		}
	}
	
	void addMmseqsAndDummies( int species ) throws IllegalValueException, CloneNotSupportedException {
		JMmseqs mmseqs = new JMmseqs(species);
		add( mmseqs );
		for( int i = 1; i < threads; i++ ) {
			add( new JSleep(mmseqs) );
		}
	}
	
	class JSleep extends FlaggedRunnable {
		JMmseqs mmseqs;
		
		JSleep( JMmseqs mmseqs ) {
			this.mmseqs=mmseqs;
		}
		
		@Override
		public void doJob() throws Exception {
			mmseqs.latch.await();
		}
		
	}
	
	/**
	 * Job for running mmseqs.
	 * 
	 * @author Jens Keilwagen
	 */
	class JMmseqs extends FlaggedRunnable {
		int speciesIndex;
		CountDownLatch latch;
		
		JMmseqs( int species ) throws IllegalValueException, CloneNotSupportedException {
			this.speciesIndex=species;
			latch = new CountDownLatch(1);
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting mmseqs for species " + speciesIndex + "\n");
			Species sp = species.get(speciesIndex);

			int exitCode=0;
			//createDB for query
			exitCode = runProcess( mmseqs, "createdb", escape(sp.cds), escape(home+speciesIndex+"/mmseqsdb"), /*"--dont-split-seq-by-len", "--dont-shuffle", "--dbtype", "1",*/ "-v", "2" );

			String tmp = home+speciesIndex+"/mmseqsdb_tmp";
			File t = new File( tmp );
			t.mkdirs();

			//search
			if( exitCode==0 ) {
				exitCode = runProcess( mmseqs, "search",
					escape(home+speciesIndex+"/mmseqsdb"),
					escape(home+"/mmseqsdb"),
					escape(home+speciesIndex+"/mmseqsdb_align.out"),
					escape(tmp),
					"-e", "100.0", "--threads", ""+threads, "-s", "8.5", "-a", "--comp-bias-corr", "0", "--max-seqs", "500", "--mask", "0", "--orf-start-mode", "1", "-v", "2" ); 
			}

			//convertalis
			if( exitCode == 0 ) {
				exitCode = runProcess( mmseqs, "convertalis",
						escape(home+speciesIndex+"/mmseqsdb"),
						escape(home+"/mmseqsdb"),
						escape(home+speciesIndex+"/mmseqsdb_align.out"),
						escape(home+speciesIndex+"/mmseqs.tabular"),
						"--threads", ""+threads, "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,empty,raw,nident,empty,empty,empty,qframe,tframe,qaln,taln,qlen,tlen", "-v", "2");
			}
			latch.countDown();
			
			//ExternalSort (and split)
			if( exitCode == 0 ) {
				File[] parts = Tools.externalSort(home+speciesIndex+"/mmseqs.tabular", 500000, threads, pipelineProtocol);
				sp.searchResults = new String[parts.length];
				for( int i = 0; i < sp.searchResults.length; i++ ) {
					sp.searchResults[i] = parts[i].getAbsolutePath();
				}
			}
			
			if( exitCode == 0 ) {
				for( int i = 0; i < sp.searchResults.length; i++ ) {
					if( !queue.isShutdown() ) {
						add(new JGeMoMa(speciesIndex, i));
					}
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
		int speciesIndex, split;
		
		JGeMoMa( int species, int split ) {
			this.speciesIndex = species;
			this.split = split;
		}
		
		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting GeMoMa split=" + split + " for species " + speciesIndex + "\n");
			
			Species sp = species.get(speciesIndex);
			ToolParameterSet params = gemomaParams.clone();
			params.getParameterForName("search results").setValue(sp.searchResults[split]);
			params.getParameterForName("cds parts").setValue(sp.cds);
			params.getParameterForName("assignment").setValue(sp.assignment);
			if( sp.protein != null ) {
				params.getParameterForName("query proteins").setValue(sp.protein);
			}
			
			SysProtocol protocol = new QuietSysProtocol();
			GeMoMa gemoma = new GeMoMa(maxSize,timeOut,maxTimeOut);
			ToolResult res = gemoma.run(params, protocol, new ProgressUpdater(), 1);
			String outDir = home + "/" + speciesIndex + "/"+split +"/";
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
			for( int sp = 0; sp <species.get(speciesIndex).searchResults.length; sp++ ) {
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
			SysProtocol protocol = new QuietSysProtocol();
			
			GeMoMaAnnotationFilter gaf = new GeMoMaAnnotationFilter();
			
			ExpandableParameterSet eps = (ExpandableParameterSet) gafParams.getParameterForName("predicted annotation").getValue();
			int i = 0;
			for( int s = 0; s < speciesCounter; s++ ) {
				Species current  = species.get(s);
				if( current.hasCDS ) {
					if( i>0 ) {
						eps.addParameterToSet();
					}
					ParameterSet ps = ((ParameterSetContainer) eps.getParameterAt(i++)).getValue();
					if( current.id != null ) {
						ps.getParameterForName("prefix").setValue(current.id);
					}
					ps.getParameterForName("gene annotation file").setValue(home+"/" + s + "/unfiltered-predictions.gff");
					if( current.weight != null ) {
						ps.getParameterForName("weight").setValue(current.weight);
					}
				}
			}
			
			ToolResult tr = gaf.run(gafParams, protocol, new ProgressUpdater(), 1);
			//res.add( tr.getRawResult()[0].getResultAt(0) );
			CLI.writeToolResults(tr, protocol, home, gaf, gafParams);
		}	
	}
	
	/**
	 * Job for running AnnotationFinalizer predicting UTRs and allowing to rename the predictions.
	 * 
	 * @author Jens Keilwagen
	 */
	class JAnnotationFinalizer extends FlaggedRunnable {

		@Override
		public void doJob() throws Exception {
			pipelineProtocol.append("starting AnnotationFinalizer\n");
			SysProtocol protocol = new QuietSysProtocol();
			
			AnnotationFinalizer af = new AnnotationFinalizer();
			
			afParams.getParameterForName("annotation").setValue(home+"/filtered_predictions.gff");

			if( rnaSeq && "YES".equals( afParams.getParameterForName("UTR").getValue() ) ) {
				setRNASeqParams(afParams, pipelineProtocol);
			}
			
			ToolResult tr = af.run(afParams, protocol, new ProgressUpdater(), 1, params, getToolVersion() );
			res.add( tr.getRawResult()[0].getResultAt(0) );
			CLI.writeToolResults(tr, protocol, home, af, gafParams);
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
			SysProtocol protocol = new QuietSysProtocol();
			ToolResult tr = extractor.run(params, protocol, new ProgressUpdater(), 1);
			ResultSet raw = tr.getRawResult()[0];
			for( int i = 2; i < raw.getNumberOfResults(); i++ ) {
				TextResult sub = (TextResult) raw.getResultAt(i);
				res.add( new TextResult( "predicted " + sub.getName(), sub.getComment(), sub.getValue(), sub.getMime(), (String) sub.getProducer(), sub.getExtendedType(), sub.getExport() ) );
			}
			CLI.writeToolResults(tr, protocol, home, extractor, params);
		}	
	}

	@Override
	public String getToolName() {
		return "GeMoMa pipeline";
	}

	@Override
	public String getShortName() {
		return "GeMoMaPipeline";
	}

	@Override
	public String getDescription() {
		return "runs the complete GeMoMa pipeline for one target species and multiple reference species";
	}

	@Override
	public String getHelpText() {
		return "**What it does**\n\nThis tool can be used to run the complete GeMoMa pipeline. The tool is multi-threaded and can utilize all compute cores on one machine, but not distributed as for instance in a compute cluster.\nIt basically runs: makeblastdb, Extract RNA-seq evidence (ERE), Extractor, FastaSplitter, tblastn, GeMoMa, cat, and GeMoMa Annotation Filter (GAF).\n\n"
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return new ResultEntry[] {
				new ResultEntry(TextResult.class, "gff", AnnotationFinalizer.defResult ),
		};
	}
	
	@Override
	public ToolResult[] getTestCases() {
		try {
			return new ToolResult[]{new ToolResult(FileManager.readFile("tests/gemoma/xml/gemomapipeline-test.xml"))};
		} catch( Exception e ) {
			e.printStackTrace();
			return null;
		}
	}
}