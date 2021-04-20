package projects.gemoma;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Locale;

import javax.script.ScriptEngine;
import javax.script.ScriptEngineManager;
import javax.script.ScriptException;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.parameters.validation.RegExpValidator;
import de.jstacs.results.ResultSet;
import de.jstacs.results.TextResult;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;
import de.jstacs.utils.DoubleList;
import projects.gemoma.GeMoMa.IntArrayComparator;

/**
 * This class allows a detailed comparison of the true and predicted annotation.
 * It computes transcript sensitivity ad transcript precision and returns a table including all true and predicted transcripts.
 * Transcripts from both annotations are compared on basis of F1 measure.
 * For each predicted transcript, best best matching true transcript is determined based on F1 measure.
 * A negative value in a F1 measure column indicates that there is a predicted transcript that matches the true transcript with
 * a F1 measure value that is the absolute value of this entry, but there is another true transcript that matches this predicted
 * transcript with an even better F1. 
 * 
 * The return table can be augmented with additional information if {@link AnnotationEvidence} is run on the true annotation before.
 * 
 * @author Jens Keilwagen
 *
 * @see AnnotationEvidence
 */
public class Analyzer extends GeMoMaModule {
	
	static int getPartsBeforePos( int pos, ArrayList<int[]> list ) {
		int anz = 0, i = 0;
		while( i < list.size() && list.get(i)[1] < pos ) {
			i++;
			anz++;
		}
		return anz;
	}

	static int getPartsAfterPos( int pos, ArrayList<int[]> list ) {
		int anz = 0, i = list.size()-1;
		while( i>=0 && pos < list.get(i)[0] ) {
			i--;
			anz++;
		}
		return anz;
	}

	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads,
			String temp) throws Exception {
		FileParameter fp = (FileParameter) parameters.getParameterForName("truth");
		String truth = fp.getValue();
		
		String feature = ((Boolean) parameters.getParameterForName("CDS").getValue()) ? "CDS" : "exon" ;
		protocol.append("selected feature: " + feature + "\n\n");
		
		SimpleParameterSet ps = (SimpleParameterSet) parameters.getParameterForName("write").getValue();
		boolean write = ps.getNumberOfParameters()>0;
		double threshold =  write ? (Double) ps.getParameterForName("common attributes").getValue() : 0;
		
		ps = (SimpleParameterSet) parameters.getParameterForName("reliable").getValue();
		String filter = Tools.prepareFilter(ps.getNumberOfParameters()==0? null : (String) ps.getParameterForName("filter").getValue() );
		ScriptEngineManager mgr = new ScriptEngineManager();
		ScriptEngine engine = mgr.getEngineByName("nashorn");
		
		String sep, eol; //TODO
		sep="\t"; eol="\n";
		//sep="\t&"; eol="\\\\\n"; //Latex
		
		//read and organize truth
		protocol.append("reading true annotation and removing duplicate transcripts\n");
		HashMap<String,int[]> attributesTruth = new HashMap<String,int[]>();
		attributesTruth.put(ALL,new int[1]);
		ArrayList<String> attTruth = new ArrayList<String>();
		HashMap<String,HashMap<String,Transcript>> res = readGFF( feature, truth, protocol, attributesTruth, attTruth );
		HashMap<String,Gene[]> genes = toGenes(res, engine, filter, protocol);

		ArrayList<String> n = new ArrayList<String>();
		ArrayList<GFFCompareStat> stats = new ArrayList<GFFCompareStat>();
		ArrayList<TextResult> tr = new ArrayList<TextResult>();
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("prediction").getValue();

		for( int i = 0; i < eps.getNumberOfParameters(); i++ ) {
			SimpleParameterSet sps = ((SimpleParameterSet)eps.getParameterAt(i).getValue());
			String predicted = (String) sps.getParameterForName("predicted annotation").getValue();
			String name = (String) sps.getParameterForName("name").getValue();
			if( name == null || name.length()==1 ) {
				name = predicted;
			}
			n.add(name);
			
			//read prediction
			protocol.append("reading predicted annotation: " + name + "\n");
			HashMap<String,int[]> attributesPrediction = new HashMap<String,int[]>();
			attributesPrediction.put(ALL,new int[1]);
			ArrayList<String> attPrediction = new ArrayList<String>();
			HashMap<String,HashMap<String,Transcript>> prediction = readGFF( feature, predicted, protocol, attributesPrediction, attPrediction );
			
			//compare
			protocol.append("comparing true and predicted annotation\n");
			HashMap<String,Transcript[]> noOverlap = compare( genes, prediction, protocol, filter.length()>0, stats );
			
			//write result
			if( write ) {
				tr.add( write( "Comparison " + i, genes, noOverlap, attributesTruth, attTruth, attributesPrediction, attPrediction, temp, threshold, filter ) );
			}
			protocol.append("\n");
		}
		
		write(protocol, filter.length()>0, n, stats, sep, eol );
		
		if( tr.size()>0 ) {
			return new ToolResult("", "", null, new ResultSet(tr), parameters, getToolName(), new Date());
		} else {
			return null;
		}
	}
	
	static final String ALL = "ALL";
	
	public TextResult write( String name, HashMap<String,Gene[]> genes, HashMap<String,Transcript[]> noOverlap, HashMap<String,int[]> attributesTruth,
			ArrayList<String> attTruth, HashMap<String,int[]> attributesPrediction,
			ArrayList<String> attPrediction, String temp, double threshold, String filter ) throws IOException {
		File out = Tools.createTempFile("Analyzer", temp);
		BufferedWriter b = new BufferedWriter( new FileWriter(out) );
		
		String[] chrom = chr.toArray(new String[0]);
		Arrays.sort(chrom);
		Gene[] both = new Gene[2];
		int[] index = new int[2];
		b.append( "chr\tstrand\tgeneId\tgene bestF1\t#transcripts\ttranscriptId");
		boolean rel = filter.length()>0;
		if( rel ) b.append( "\treliable (" + filter + ")" );
		b.append("\ttranscript start position\ttranscript end position" );
		double all = attributesTruth.get(ALL)[0];
		int anzTruth=0;
		for( int i = 0; i < attTruth.size(); i++ ) {
			if( attributesTruth.get(attTruth.get(i))[1]/all >= threshold ) {
				b.append( "\ttranscript " + attTruth.get(i));
				anzTruth++;
			}
		}
		b.append( "\ttranscript bestF1\ttranscript F1"
				
				//stats
				+ "\tstart\tend\tfeature difference\tfirst feature\tlast feature\ttp\tfn\tfp"
				+ "\tadditional upstream features truth\tadditional internal features truth\tadditional downstream features truth"
				+ "\tadditional upstream features prediction\tadditional internal features prediction\tadditional downstream features prediction"
				+ "\tintron retention in truth\tintron retention in prediction"
				+ "\tperfect features"
				
				+ "\tpredictionId\tprediction start position\tprediction end position" );
		all = attributesPrediction.get(ALL)[0];
		for( int i = 0; i < attPrediction.size(); i++ ) {
			if( attributesPrediction.get(attPrediction.get(i))[1]/all >= threshold ) {
				b.append( "\tprediction " + attPrediction.get(i));
			}
		}
		b.newLine();
		for( int i = 0; i < chrom.length; i++ ) {
			Gene[] pos = genes.get(chrom[i]+"+");
			Gene[] neg = genes.get(chrom[i]+"-");
			Transcript[] noOv = noOverlap.get(chrom[i]);
			
			Arrays.fill(index, 0);
			int pos_max = pos==null?0:pos.length, neg_max= neg==null?0:neg.length,
					noOv_j = 0, noOv_max= noOv==null?0:noOv.length;
			while( index[0] < pos_max || index[1] < neg_max || noOv_j < noOv_max ) {
				both[0] = index[0] < pos_max ? pos[index[0]] : null;
				both[1] = index[1] < neg_max ? neg[index[1]] : null;
				
				Transcript o = noOv_j < noOv_max ? noOv[noOv_j] : null;

				int best;
				if( both[0]==null && both[1]==null ) {
					best = -1;
				} else {
					if( (both[0]!=null && both[1]!= null && both[0].compareTo(both[1])< 0) || both[1]==null ) {
						best=0;
					} else {
						best=1;
					}
				}

//System.out.println(best + "\t" + index[0] + "/" + pos_max + "\t" + index[1] + "/" + neg_max + "\t" + noOv_j + "/" +noOv_max );				
				if( o==null || (best >= 0 && o != null && both[best].min < o.min) ) {
					both[best].write( attributesTruth, attTruth, attributesPrediction, attPrediction, b, threshold, rel );
					index[best]++;
				} else {
					
					noOv_j++;
					tBuff.delete(0,tBuff.length());
					tBuff.append( o.chr + "\t" + o.strand + "\t\tNA\t\t\t\t" + (rel?"\t":"") );
					for( int j = 0; j < anzTruth; j++ ) {
						tBuff.append( "\t" );
					}
					tBuff.append( "\tNA\tNA" );
					//last part is for additional fields: start, end, number of exons, ...
					for( int k = 0; k < 17; k++ ) {
						tBuff.append( "\t" );
					}
					o.appendAttributes( tBuff, attributesPrediction, attPrediction, threshold, false ); 
					b.append( tBuff.toString() );
					b.newLine();
				}
			}
		}
		b.close();
		
		return new TextResult(name, "Result", new FileParameter.FileRepresentation(out.getAbsolutePath()), "tabular", getToolName(), null, true);
	}
	
	public HashMap<String,Transcript[]> compare( HashMap<String,Gene[]> genes, HashMap<String,HashMap<String,Transcript>> prediction, Protocol protocol, boolean rel, ArrayList<GFFCompareStat> stats ) {
		HashMap<String,Transcript[]> noOverlap = new HashMap<String,Transcript[]>();
		ArrayList<Transcript> current = new ArrayList<Transcript>();
		Iterator<String> it = chr.iterator();
		BitSet predictedRegion = new BitSet(), trueRegion = new BitSet(), help = new BitSet();
		GFFCompareStat stat = new GFFCompareStat();
		HashMap<String,Double> predGenes = new HashMap<String,Double>();
		
		while( it.hasNext() ) {
			String chrom = it.next();
			protocol.append( chrom + "\n" );
			current.clear();
			for( int i = 0; i < 2; i++ ) {
				String c = chrom+(i==0?"+":"-");
				Gene[] g = genes.get(c);
				HashMap<String,Transcript> pred = prediction.get(c);
				
				if( g != null ) {
					for( int j=0; j < g.length; j++ ) {
						stat.anzTruth += g[j].t.size();
						boolean r = false;
						for( int k = 0; k < g[j].t.size(); k++ ) {
							Transcript t = g[j].t.get(k);
							//delete old
							t.clear();
							
							if( t.reliable ) {
								stat.anzRTruth++;
								r=true;
							}
						}
						if( r ) stat.anzRTruthG++;
					}
					stat.anzTruthG +=g.length;
				}
				if( pred!= null ) {
					stat.anzPrediction += pred.size();
					Iterator<Transcript> iter = pred.values().iterator();
					while( iter.hasNext() ) {
						predGenes.put( iter.next().parent, 0d );
					}
					if( g!= null ) {
						Transcript[] predictions = pred.values().toArray(new Transcript[0]);
						Arrays.sort(predictions);

						int[] cumGEnd = new int[g.length];
						int last = -1000;
						for( int j = 0; j < g.length; j++ ) {
							last = cumGEnd[j] = Math.max(last, g[j].max);
						}
						int gIndex=0;
						for( int j = 0; j < predictions.length; j++ ) {
							predictions[j].sortParts();
							//skip genes that end before the current prediction starts
							while( gIndex < cumGEnd.length && cumGEnd[gIndex] < predictions[j].min ) {
								gIndex++;
							}
							//check all potential genes the could overlap with the prediction
							double bestF1 = -1;
							int bestG = -1, bestT = -1;
							int idx = gIndex;
							while( idx < g.length && g[idx].min < predictions[j].max ) {
								for( int k = 0; k < g[idx].t.size(); k++ ) {
									Transcript currentTranscript = g[idx].t.get(k);
									int start = Math.min( predictions[j].min, currentTranscript.min )-1;
									//int end = Math.max( predictions[j].max, currentTranscript.max )+1;
									predictions[j].set( predictedRegion, start );
									currentTranscript.set( trueRegion, start );
									
									double tp = count( help, predictedRegion, trueRegion, false); //tp
									double fn = count( help, predictedRegion, trueRegion, true); //fn
									double fp = count( help, trueRegion, predictedRegion, true); //fp
			
									double f1 = 2d*tp/(double)(2d*tp+fn+fp);
									if( f1 > 0 ) {
										currentTranscript.setOverlap(f1);
									}
									if( f1 > bestF1 ) {
										bestF1=f1;
										bestG=idx;
										bestT=k;
									}
								}
								idx++;
							}
							
							if( predGenes.get(predictions[j].parent) < bestF1 ) {
								predGenes.put(predictions[j].parent, bestF1);
							}
							
							if( bestF1<=0 ) {
								current.add(predictions[j]);
							} else {
								int start = Math.min( predictions[j].min, g[bestG].t.get(bestT).min )-1;
								
								predictions[j].set( predictedRegion, start );
								g[bestG].t.get(bestT).set( trueRegion, start );
								
								double tp = count( help, predictedRegion, trueRegion, false); //tp
								double fn = count( help, predictedRegion, trueRegion, true); //fn
								double fp = count( help, trueRegion, predictedRegion, true); //fp
								
								g[bestG].t.get(bestT).set( predictions[j], bestF1, tp, fn, fp );
								
								if( bestF1 == 1 ) {
									stat.anzPerfect++;
									if( g[bestG].t.get(bestT).reliable ) stat.anzRPerfect++;
								}
							}
						}
						
						for( int j=0; j < g.length; j++ ) {
							if( g[j].getBestF1()==1 ) {
								stat.anzPerfectG++;
							}
							if( rel ) {
								boolean r = false;
								for( int k = 0; k < g[j].t.size(); k++ ) {
									Transcript t = g[j].t.get(k);
									if( t.reliable && t.bestF1==1 ) {
										r=true;
									}
								}
								if( r ) stat.anzRPerfectG++;
							}
						}
					} else {
						current.addAll(pred.values());
					}
				}
			}

			Collections.sort(current);
			noOverlap.put(chrom, current.toArray(new Transcript[0]) );
		}

		stat.anzPredG=predGenes.size();
		Iterator<Double> f = predGenes.values().iterator();
		while( f.hasNext() ) {
			if( f.next() == 1d ) {
				stat.anzPerfectG2++;
			}
		}
		stats.add(stat);
		
		return noOverlap;
	}
	
	class GFFCompareStat {
		int anzTruth = 0, anzRTruth=0, anzPrediction=0, anzTruthG=0, anzRTruthG=0, anzPredG;
		double anzPerfect=0, anzRPerfect=0, anzPerfectG=0, anzPerfectG2=0, anzRPerfectG=0;		
	}
	
	static void write( Protocol protocol, boolean rel, ArrayList<String> name, ArrayList<GFFCompareStat> list, String sep, String eol ) {
		Locale l = Locale.US;
		int prec=2;
		
		NumberFormat nf1 = NumberFormat.getInstance(l);
		nf1.setMinimumFractionDigits(prec);
		nf1.setMaximumFractionDigits(prec);
		NumberFormat nf2 = NumberFormat.getInstance(l);
		nf2.setMaximumFractionDigits(0);
		
		protocol.append( "number of true transcripts");
		int val = list.get(0).anzTruth;
		for( GFFCompareStat stat: list ) if(val!=stat.anzTruth) throw new IllegalArgumentException("anzTruth corrupted") ;
		protocol.append( sep + nf2.format(val) + eol );
		protocol.append( "number of true genes" );
		int val2 = list.get(0).anzTruthG;
		for( GFFCompareStat stat: list ) if(val2!=stat.anzTruthG) throw new IllegalArgumentException("anzTruthG corrupted") ;
		protocol.append( sep + nf2.format(val2) + eol );
		protocol.append( "number of true transcripts per true gene" + sep + nf1.format(val/(double)val2) + eol );
		
		if( rel ) {
			protocol.append( "number of reliable true transcripts" );
			val = list.get(0).anzRTruth;
			for( GFFCompareStat stat: list ) if(val!=stat.anzRTruth) throw new IllegalArgumentException("anzRTruthG corrupted") ;
			protocol.append( sep + nf2.format(val)  + " (" + nf1.format(val/(double)list.get(0).anzTruth*100) + "%)" + eol );
			protocol.append( "number of reliable true genes" );
			val = list.get(0).anzRTruthG;
			for( GFFCompareStat stat: list ) if(val!=stat.anzRTruthG) throw new IllegalArgumentException("anzRTruthG corrupted") ;
			protocol.append( sep + nf2.format(val) + " (" + nf1.format(val/(double)list.get(0).anzTruthG*100) + "%)"  + eol );
		}
		
		protocol.append( eol+eol+"category");
		for( String n: name ) protocol.append( sep + n );
		protocol.append( eol+eol );
		
		protocol.append( "number of predicted transcripts" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzPrediction) );
		protocol.append( eol );
		protocol.append( "number of predicted genes" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzPredG) );
		protocol.append( eol );
		protocol.append( "number of predicted transcripts per predicted gene" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzPrediction/(double)stat.anzPredG) );
		protocol.append( eol );
		
		protocol.append( eol + "number of perfectly predicted transcripts" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzPerfect) );
		protocol.append( eol );
		if( rel ) {
			protocol.append( "number of perfectly predicted reliable true transcripts" );
			for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzRPerfect) );
			protocol.append( eol );
		}

		protocol.append( eol + "number of true genes with at least one perfectly predicted transcript" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzPerfectG) );
		protocol.append( eol );
		protocol.append( "number of predicted genes with at least one perfectly predicted transcript" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzPerfectG2) );
		protocol.append( eol );
		if( rel ) {
			protocol.append( "number of reliable true genes with at least one perfectly predicted transcript" );
			for( GFFCompareStat stat: list ) protocol.append( sep + nf2.format(stat.anzRPerfectG) );
			protocol.append( eol );
		}
		
		protocol.append( eol + "transcript sensitivity");
		for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzPerfect / stat.anzTruth*100) );
		protocol.append( eol );
		protocol.append( "transcript precision" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzPerfect / stat.anzPrediction*100) );
		protocol.append(eol );
		if( rel ) {
			protocol.append( "reliable transcript sensitivity" );
			for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzRPerfect / stat.anzRTruth*100) );
			protocol.append( eol );
		}
		
		protocol.append( eol + "gene sensitivity" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzPerfectG / stat.anzTruthG*100) );
		protocol.append( eol );
		protocol.append( "gene precision" );
		for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzPerfectG2 / stat.anzPredG*100) );
		protocol.append( eol );
		if( rel ) {
			protocol.append( "reliable gene sensitivity" );
			for( GFFCompareStat stat: list ) protocol.append( sep + nf1.format(stat.anzRPerfectG / stat.anzRTruthG *100) );
			protocol.append( eol );
		}
	}

	
	private static int count( BitSet help, BitSet b1, BitSet b2, boolean flipFirst ) {
		help.clear();
		help.or(b1);
		if( flipFirst ) {
			help.flip(0,help.size());
		}
		help.and(b2);
		return help.cardinality();
	}
	
	public HashMap<String,Gene[]> toGenes( HashMap<String,HashMap<String,Transcript>> res, ScriptEngine engine, String filter, Protocol protocol ) throws ScriptException {
		Iterator<String> it = res.keySet().iterator();
		HashMap<String,Gene> g = new HashMap<String,Gene>();
		HashMap<String,Gene[]> sortedGenes = new HashMap<String,Gene[]>();
		int duplicate=0;
		while( it.hasNext() ) {
			String key = it.next();
			HashMap<String,Transcript> current = res.get(key);
			//System.out.println(key + "\t" + current.size() );
			
			//to genes
			g.clear();
			Iterator<String> tit = current.keySet().iterator();
			while( tit.hasNext() ) {
				String tkey = tit.next();
				Transcript t = current.get(tkey);
				if( t.parts.size()>0 ) {
					if( filter.length()>0 ) t.reliable = Tools.filter(engine, filter, t.hash);
					t.sortParts();
					if( t.parent == null ) {
						t.parent = t.id+".gene";
					}
					Gene gene = g.get(t.parent);
					if( gene == null ) {
						gene = new Gene(t.parent);
						g.put(t.parent, gene);
					}
					gene.addTranscript(t);
				}
			}
			Gene[] gArray=g.values().toArray(new Gene[0]);
			Arrays.sort(gArray);
			for( int i = 0; i < gArray.length; i++ ) {
				duplicate+=gArray[i].eliminateDuplicates();
				Collections.sort( gArray[i].t, Analyzer.TranscriptNameComparator.SINGLETON );
			}
			
			sortedGenes.put( key, gArray );
			//System.out.println(key + "\t" + gArray.length );
		}
		if( duplicate>0 ) protocol.append("removed " + duplicate + " duplicates\n\n" );
		//System.out.println(key + "\t" + gArray.length );
		return sortedGenes;
	}

	public HashMap<String,HashMap<String,Transcript>> readGFF( String feature, String fName, Protocol protocol, HashMap<String,int[]> attributes, ArrayList<String> att ) throws IOException {
		BufferedReader r = new BufferedReader( Tools.openGzOrPlain(fName) );
		HashMap<String,HashMap<String,Transcript>> res = new HashMap<String,HashMap<String,Transcript>>();
		HashMap<String,Transcript> current = new HashMap<String,Transcript>();
		String line;
		String[] split;
		while( (line=r.readLine()) != null ) {
			if( line.equalsIgnoreCase("##FASTA") ) {//http://gmod.org/wiki/GFF3#GFF3_Sequence_Section 
				protocol.append("Stop reading the annotation file because of '##FASTA'\n"); 
				break;  
			}
			if( line.length() == 0 || line.startsWith("#") ) continue;

			split = line.split("\t");
			if( split.length!=9 ) {
				throw new IllegalArgumentException("This line does not seem to be a valid (tab-delimited) line of the annotation with 9 columns: " + line);
			}
			
			boolean relevant = split[2].equals(feature);
			if( relevant || split[2].equals("mRNA") || split[2].equals("transcript") ) {
				chr.add(split[0]);
				current = res.get(split[0]+split[6]);
				if( current == null ) {
					current = new HashMap<String,Transcript>();
					res.put(split[0]+split[6], current);
				}
				
				int idx1, idx2;
				if( relevant ) {
					idx1 = split[8].indexOf("Parent=")+7;
					idx2 = split[8].indexOf(';',idx1);
					if( idx2<0 ) idx2=split[8].length();
				} else {
					idx1 = split[8].indexOf("ID=")+3;
					idx2 = split[8].indexOf(';',idx1);
					if( idx2<0 ) idx2=split[8].length();
				}
				if( idx1<0 || idx2<0 ) {
					throw new IllegalArgumentException(line);
				}
				String tID = split[8].substring(idx1,idx2);
				Transcript t = current.get( tID );
				if( t == null ) {
					t = new Transcript(split[0], split[6], tID );
					current.put( tID, t );
				}
				if( relevant ) {
					t.addFeature(split[3], split[4]);
				} else {
					t.addInfo(split[8],attributes,att);
				}
			}
		}
		r.close();
		
		int empty=0;
		Iterator<String> it = res.keySet().iterator();
		while( it.hasNext() ) {
			current = res.get(it.next());
			String[] ids = current.keySet().toArray(new String[0]);
			for( int i = 0; i < ids.length; i++ ) {
				Transcript t = current.get(ids[i]);
				if( t.parts.size() == 0 ) {
					empty++;
					current.remove(ids[i]);
				}
			}
		}
		if( empty>0 ) protocol.append("removed " + empty + " empty transcripts\n" );
		return res;
	}
	
	HashSet<String> chr = new HashSet<String>();
	StringBuffer tBuff = new StringBuffer();
	
	/**
	 * This class represents a transcript - either a true or a predicted transcript.
	 * 
	 * @author Jens Keilwagen
	 */
	class Transcript implements Comparable<Transcript> {
		String chr, id, parent;
		char strand;
		HashMap<String,String> hash;
		ArrayList<int[]> parts;
		int min, max;
		
		ArrayList<Transcript> overlap;
		DoubleList overlapF1;
		ArrayList<String> overlapInfo;
		double bestF1;
		
		boolean reliable;
		
		Transcript( String chr, String strand, String id ) {
			this.chr=chr;
			this.strand=strand.charAt(0);
			this.id = id;
			parent = null;
			hash=null;
			parts = new ArrayList<int[]>();
			min = Integer.MAX_VALUE;
			max = -100;
			clear();
			reliable=true;
		}

		public void sortParts() {
			Collections.sort( parts, IntArrayComparator.comparator[2] );
		}

		public void set( BitSet region, int start ) {
			region.clear();
			for( int i = 0; i < parts.size(); i++ ) {
				int[] interval = parts.get(i);
				region.set(interval[0]-start, interval[1]-start+1);
			}
		}
		
		public void clear() {
			bestF1=0;
			if( overlap!=null ) overlap.clear();
			if( overlapF1!=null ) overlapF1.clear();
			if( overlapInfo!=null ) overlapInfo.clear();
		}
		
		public void set(Transcript transcript, double f1, double tp, double fn, double fp) {
			if( overlap == null ) {
				overlap = new ArrayList<Transcript>();
				overlapF1 = new DoubleList();
				overlapInfo = new ArrayList<String>();
			}
			overlap.add(transcript);
			overlapF1.add(f1);
			String s;
			
			//min max
			if( strand=='+' ) {
				s = (min==transcript.min)  + "\t" + (max==transcript.max);
			} else {
				s = (max==transcript.max)  + "\t" + (min==transcript.min);
			}
			//difference
			s+= "\t" + (parts.size() - transcript.parts.size());
			if( strand=='+' ) {
				s += "\t" + (IntArrayComparator.comparator[0].compare(parts.get(0), transcript.parts.get(0))==0)
						+ "\t" + (IntArrayComparator.comparator[0].compare(parts.get(parts.size()-1), transcript.parts.get(transcript.parts.size()-1))==0);
			} else {
				s += "\t" + (IntArrayComparator.comparator[0].compare(parts.get(parts.size()-1), transcript.parts.get(transcript.parts.size()-1))==0)
						+ "\t" + (IntArrayComparator.comparator[0].compare(parts.get(0), transcript.parts.get(0))==0);
			}
			
			//tn, fn, fp
			s+= "\t" + tp + "\t" + fn + "\t" + fp;
			
			//upstream/downstream feature
			int a = getPartsBeforePos( transcript.parts.get(0)[0], parts );
			int b = getPartsAfterPos( transcript.parts.get(transcript.parts.size()-1)[1], parts );
			int aa = getPartsBeforePos( parts.get(0)[0], transcript.parts );
			int bb = getPartsAfterPos( parts.get(parts.size()-1)[1], transcript.parts );
			
			//intron retention
			int i = 0, j = 0;
			int intronRetentionTruth = 0, intronRetentionPrediction = 0;
			int addIntFeatureTruth=0, addIntFeaturePrediction=0;
			while( i < parts.size() && j < transcript.parts.size() ) {
				int[] tru = parts.get(i);
				int[] pred = transcript.parts.get(j);
				if( tru[1]<pred[0] ) { //tru lies before pred
					if( j>0 && transcript.parts.get(j-1)[1]<tru[0] ) {
						addIntFeatureTruth++;
					}
					i++; 
				} else if( pred[1]<tru[0] ) { //pred lies before tru
					if( i>0 && parts.get(i-1)[1]<pred[0] ) {
						addIntFeaturePrediction++;
					}
					j++;
				} else { //overlap between t and p
					if( pred[1]<tru[1] ) { //truth is longer
						if( j+1 < transcript.parts.size() && transcript.parts.get(j+1)[0] < tru[1] ) {
							intronRetentionTruth++;
						}
						j++;
					} else if ( tru[1]<pred[1] ) { //prediction is longer
						if( i+1 < parts.size() && parts.get(i+1)[0] < pred[1] ) {
							intronRetentionPrediction++;
						}						
						i++;
					} else { //t[1]==p[1]
						i++;
						j++;
					}
				}
			}
			if( strand=='+' ) {
				s += "\t" + a + "\t"+addIntFeatureTruth+ "\t" + b + "\t" + aa +"\t"+addIntFeaturePrediction+ "\t" + bb;
			} else {
				s += "\t" + b + "\t"+addIntFeatureTruth+ "\t" + a + "\t" + bb +"\t"+addIntFeaturePrediction+ "\t" + aa;
			}

			s += "\t"+intronRetentionTruth+"\t"+intronRetentionPrediction;
			
			//perfect exons
			i = j = 0;
			int anz=0;
			while( i < parts.size() && j < transcript.parts.size() ) {
				int comp = IntArrayComparator.comparator[0].compare(parts.get(i), transcript.parts.get(j));
				if( comp == 0 ) {
					i++;
					j++;
					anz++;
				} else if( comp < 0 ) {
					i++;
				} else {
					j++;
				}
			}
			s += "\t" + anz;

			
			//TODO extend?
			
			overlapInfo.add(s);
		}
		
		public void setOverlap(double f1) {
			if( bestF1 < f1 ) {
				bestF1=f1;
			}
		}
		
		public double bestF1() {
			if( overlap == null || overlap.size()==0 ) {
				if( bestF1>0 ) {
					return -bestF1;
				} else {
					return 0;
				}
			} else {
				return overlapF1.max(0, overlapF1.length());
			}
		}

		public void addInfo( String attributes, HashMap<String,int[]> attributesHash, ArrayList<String> att ) {
			if( hash!= null || parent != null ) {
				throw new IllegalArgumentException("Could not set attributes ("+hash+") or parent ("+parent+"): " + attributes );
			}
			attributesHash.get(ALL)[0]++;
			String[] split = attributes.split(";");
			hash = new HashMap<String, String>();
			for( String s: split ) {
				s = s.trim();
				if( s.length()== 0 ) continue;
				int index = s.indexOf('=');
				if( index < 0 ) {
					throw new IllegalArgumentException("Illegal attributes: " + s);
				}
				String key = s.substring(0,index);
				String value = s.substring(index+1);
				
				if( key.equals("ID") || key.charAt(0)=='#' ) {
					//ignore
				} else if( key.equals("Parent") ) {
					parent=value;
				} else {
					int[] h = attributesHash.get(key);
					if( h == null ) {
						h=new int[] {att.size(),0};
						attributesHash.put(key,h);
						att.add(key);
					}
					h[1]++;

					if( value.length()>100 ) value="<TOO_LONG>";
					hash.put(key, value);
				}
			}
		}
		
		public void addFeature( String start, String end ) {
			int s = Integer.parseInt(start);
			int e = Integer.parseInt(end);
			parts.add(new int[] {s,e} );
			min=Math.min(min, s);
			max=Math.max(max, e);
		}
		
		public void write( StringBuffer res, String geneID, double bestF1, int size, HashMap<String,int[]> attributesTruth, ArrayList<String> attTruth, HashMap<String,int[]> attributesPredcition, ArrayList<String> attPredcition, double threshold, boolean rel ) {
			res.delete(0, res.length());
			res.append( chr + "\t" + strand + "\t" + geneID + "\t" + bestF1 + "\t" + size );
			appendAttributes(res, attributesTruth, attTruth, threshold, rel);
			res.append( "\t" + bestF1() );
			String prefix = res.toString();
			if( overlap!= null && overlap.size()>0 ) {
				for( int i = 0; i < overlap.size(); i++ ) {
					if( i>0 ) res.append( prefix );
					res.append("\t" + overlapF1.get(i) + "\t" + overlapInfo.get(i));
					overlap.get(i).appendAttributes(res, attributesPredcition, attPredcition, threshold, false );
					res.append("\n");
				}
			} else {
				res.append( "\tNA\n");
			}
		}
		
		public void appendAttributes( StringBuffer res, HashMap<String,int[]> attributes, ArrayList<String> att, double threshold, boolean rel ) {
			res.append( "\t" + id + (rel ? "\t"+reliable : "") + "\t" + min + "\t" + max );
			double all = attributes.get(ALL)[0];
			for( int i = 0; i < att.size(); i++ ) {
				String key = att.get(i);
				if( attributes.get(key)[1]/all>threshold ) {
					String info = hash== null ? null : hash.get(key);
					res.append( "\t" + (info==null ? "" : info) );
				}
			}
		}

		@Override
		public int compareTo(Transcript t) {
			int d = Integer.compare(min, t.min);
			if( d == 0 ) {
				d = Integer.compare(max, t.max);
				if( d == 0 ) {
					d = id.compareTo(t.id);
				}
			}
			return d;
		}
	}
	
	/**
	 * This class represents a gene which is a collection of true transcripts.
	 * 
	 * @author Jens Keilwagen
	 */
	class Gene implements Comparable<Gene> {
		String id;
		int min, max;
		ArrayList<Transcript> t;
		
		Gene( String id ) {
			this.id=id;
			t = new ArrayList<Transcript>();
			min = Integer.MAX_VALUE;
			max = -100;
		}
		
		public int eliminateDuplicates() {
			Collections.sort( t, TranscriptComparator.SINGLETON );
			int dupl=0;
			for( int i = t.size()-1; i>0; i--) {
				Transcript current = t.get(i);
				Transcript next = t.get(i-1);
				if( TranscriptComparator.SINGLETON.compare(current, next) == 0 ) {
					//TODO add alternative names?
					if( TranscriptNameComparator.SINGLETON.compare(current, next) > 0 ) {
						t.remove(i);
						next.id += ","+current.id;
					} else {
						t.remove(i-1);
						current.id += ","+next.id;
					}
					dupl++;
				}
			}
			return dupl;
		}

		void addTranscript( Transcript transcript ) {
			if( !transcript.parent.equals(id) ) throw new IllegalArgumentException( id + "\t" + transcript.parent );
			t.add(transcript);
			min = Math.min(min, transcript.min);
			max = Math.max(max, transcript.max);
		}

		@Override
		public int compareTo(Gene g) {
			return Integer.compare(min, g.min);
		}
		
		public double getBestF1() {
			double bestF1 = -2;
			for( int i = 0; i < t.size(); i++ ) {
				bestF1 = Math.max( bestF1,  t.get(i).bestF1() );
			}
			return bestF1;
		}
		
		public void write( HashMap<String,int[]> attributesTruth, ArrayList<String> attTruth,
				HashMap<String,int[]> attributesPrediction, ArrayList<String> attPrediction, BufferedWriter b, double threshold, boolean rel ) throws IOException {
			double bestF1 = getBestF1();
			for( int i = 0; i < t.size(); i++ ) {
				Transcript trans = t.get(i);
				trans.write(tBuff, id, bestF1, t.size(), attributesTruth, attTruth, attributesPrediction, attPrediction, threshold, rel);
				b.append( tBuff.toString() );
			}
		}
	}
	
	/**
	 * This {@link Comparator} allows to compare transcripts by name (=id).
	 * 
	 * @author Jens Keilwagen
	 */
	static class TranscriptNameComparator implements Comparator<Transcript> {

		public final static TranscriptNameComparator SINGLETON = new TranscriptNameComparator();
		
		@Override
		public int compare(Transcript a, Transcript b) {
			return a.id.compareTo(b.id);
		}
	}
	
	/**
	 * This {@link Comparator} allows to compare transcripts by CDS.
	 * 
	 * @author Jens Keilwagen
	 */
	static class TranscriptComparator implements Comparator<Transcript> {
		
		public final static TranscriptComparator SINGLETON = new TranscriptComparator();
		
		@Override
		public int compare(Transcript a, Transcript b) {
			int diff = a.parts.size() - b.parts.size();
			if( diff == 0 ) {
				int i = 0;
				while( i < a.parts.size() && (diff=IntArrayComparator.comparator[0].compare(a.parts.get(i),b.parts.get(i))) == 0 ) {
					i++;
				}
			}
			return diff;
		}
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getToolName(),
				new FileParameter("truth", "the true annotation", "gff,gff3,gff.gz,gff3.gz", true, new FileExistsValidator()),
				new ParameterSetContainer( "prediction", "", new ExpandableParameterSet( new SimpleParameterSet(		
						new SimpleParameter(DataType.STRING,"name", "can be used to distinguish different predictions", false, new RegExpValidator("\\w*") ),
						new FileParameter( "predicted annotation", "GFF file containing the predicted annotation", "gff,gff3,gff.gz,gff3.gz",  true, new FileExistsValidator(), true )
					), "gene annotations", "", 1 ) ),
				//new FileParameter("prediction", "the predicted annotation", "gff,gff3", true, new FileExistsValidator()),
				new SimpleParameter( DataType.BOOLEAN, "CDS", "if true CDS features are used otherwise exon features", true, true),
				new SelectionParameter( DataType.PARAMETERSET, 
						new String[]{"NO","YES"},
						new Object[]{
							new SimpleParameterSet(),
							new SimpleParameterSet(	new SimpleParameter( DataType.DOUBLE, "common attributes", "Only gff attributes of mRNAs are included in the result table, that can be found in the given portion of all mRNAs. Attributes and their portion are handled independently for truth and prediction. This parameter allows to choose between a more informative table or compact table.", true, new NumberValidator<Double>(0d, 1d), 0.5 ) )
						},
						"write", "write detailed table comparing the true and the predicted annotation", true ),
				new SelectionParameter( DataType.PARAMETERSET, 
						new String[]{"NO","YES"},
						new Object[]{
							new SimpleParameterSet(),
							new SimpleParameterSet(	
								new SimpleParameter(DataType.STRING,"filter","A filter for deciding which transcript from the truth are reliable or not. The filter is applied to the GFF attributes of the truth. You probably need to run AnnotationEvidence on the truth GFF. "
									+ "The default filter decides based on the completeness of the prediction (start=='M' and stop=='*'), no premature stop codons (nps==0), RNA-seq coverage (tpc==1) and intron evidence (isNaN(tie) or tie==1).",
									false, new RegExpValidator("[a-zA-Z 0-9\\.()><=!'\\-\\+\\*\\/]*"), "start=='M' and stop=='*' and nps==0 and (tpc==1 and (isNaN(tie) or tie==1))"
								)
							)
						},
						"reliable", "additionally evaluate sensitivity for reliable transcripts", true
				)
			);
		} catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return getShortName();
	}

	@Override
	public String getShortName() {
		return getClass().getSimpleName();
	}

	@Override
	public String getDescription() {
		return "compares true and predicted annotation";
	}

	@Override
	public String getHelpText() {
		return "This tools allows to compare true annotation with predicted annotation as it is frequently done in benchmark studies."
				+ " Furthermore, it can return a detailed table comparing true annotation and predicted annotation which might help to identify systematical errors or biases in the predictions."
				+ " Hence, this tool might help to detect weaknesses of the prediction algorithm.\n\n"
				+ "True and predicted transcripts are evaluated based on nucleotide F1 measure."
				+ " For each predicted transcript, the true transcript with highest nucleotide F1 measure is listed."
				+ " Also true and predicted transcripts are listed that do not overlap with any transcript from the predicted and true annotation, respectively."
				+ " The table contains the attributes of the true and the predicted annotation besides some additional columns allowing to easily filter interesting examples and to do statistics.\n\n"
				+ "The evaluation can be based on CDS (default) or exon features."
				+ " The tool also reports sensitivity and precision for the categories gene and transcript."
			+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}
}