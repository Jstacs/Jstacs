package projects.gemoma;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.FileExistsValidator;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * The {@link BUSCORecomputer} allows to compute BUSCO statistics based on genes rather than on transcripts.
 * Based on the full BUSCO table for transcripts and a second table that has one column for gene ID and a second column for transcript ID, the BUSCO stistics for genes are computed.
 * 
 * @author Jens Keilwagen
 * 
 * @see Extractor
 */
public class BUSCORecomputer extends GeMoMaModule {
	
	public static String rem = "<REMAINING>";
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		ExpandableParameterSet eps = (ExpandableParameterSet) parameters.getParameterForName("subgenomes").getValue();
		int poly=eps.getNumberOfParameters();
		String[] regex;
		if( poly == 0 ) {
			regex = new String[0];
			poly = 1;
		} else {
			regex=new String[poly];
			for( int i = 0;i < poly; i++ ) {
				regex[i] = (String) ((SimpleParameterSet) eps.getParameterAt(i).getValue()).getParameterAt(0).getValue();
			}
		}
		
		FileParameter fp = (FileParameter) parameters.getParameterForName("BUSCO");
		String busco = fp.getValue();
		
		fp = (FileParameter) parameters.getParameterForName("IDs");
		String genTranscript = fp.getValue();
		
		HashMap<String,String> trans2gene = new HashMap<String,String>();
		HashMap<String,Integer> gene2subgenome = new HashMap<String,Integer>();
		BufferedReader r = new BufferedReader( new FileReader( genTranscript ) );
		int[] num = new int[poly+1];
		String line;
		while( (line=r.readLine()).charAt(0)=='#' );
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			trans2gene.put(split[1], split[0]);
			if( !gene2subgenome.containsKey(split[0]) ) {
				int match=-1;
				for( int j = 0; j < regex.length; j++ ) {
					if( split[4].matches(regex[j]) ) {
						if( match<0 ) {
							match=j;
						} else {
							throw new IllegalArgumentException(split[4] +" matches multiple regular expressions: " + regex[match] + " and " + regex[j] );
						}
					}
				}
				if( match<0 ) {
					match = regex.length;
					if( poly==regex.length ) {
						poly++;
					}
				}
				gene2subgenome.put(split[0], match);
				num[match]++;
			} else {
				//TODO?
			}
		}
		r.close();
		protocol.append("subgenome\t#transcripts\n");
		for( int c = 0; c < poly; c++ ) {
			protocol.append( (c < regex.length ? regex[c] : (regex.length==0?"":rem)) + "\t" + num[c] +"\n");
		}
		protocol.append("\n");
		
		int[][] stat = new int[poly][4];
		double all=0;
		r = new BufferedReader( new FileReader( busco ) );
		String old = null;
		HashSet<String>[] hash = new HashSet[poly];
		for( int c = 0; c < poly; c++ ) {
			hash[c] = new HashSet<String>();
		}
		int anz=0;
		while( (line=r.readLine()).charAt(0)=='#' ) {
			if( anz < 3 ) protocol.append(line+"\n");
			anz++;
		}
		protocol.append( "# " + "BUSCORecomputer\n" );
		do {
			String[] split = line.split("\t");
			
			if( old!= null && !split[0].equals(old) ) {
				add(stat,hash);
				all++;
				old=null;
			}
			
			int v = getIndex(split[1]);
			if( v==3 ) {
				/*for( int c = 0; c < poly; c++ ) {
					stat[c][v]++;
				}*/
				all++;
			} else {
				String gene = trans2gene.get(split[2]);
				if( gene == null ) {
					gene = split[2];
					protocol.append("Warning no gene found for transcript: " + gene + "\n");
				}
				Integer sub = gene2subgenome.get(gene);
				if( v>= 0 ) {
					stat[sub][v]++;
					/*
					for( int c = 0; c < poly; c++ ) {
						if( c==sub ) {
							stat[sub][v]++;
						} else {
							stat[c][3]++;
						}
					}*/
					all++;
				} else {
					if( old == null ) {
						old=split[0];
					}
					hash[sub].add(gene);
				}
			}
		} while( (line=r.readLine()) != null );
		if( old!= null ) {
			add(stat,hash);
			all++;
		}
		r.close();
		
		NumberFormat nf = NumberFormat.getInstance(Locale.US);
		nf.setMaximumFractionDigits(1);
		int a = (int) all;
		for( int c = 0; c < poly; c++ ) {
			stat[c][3] = (int)all - (stat[c][0]+stat[c][1]+stat[c][2]);
			
			protocol.append( c < regex.length ? (regex[c]+"\t") : (regex.length==0?"":(rem+"\t")));
			protocol.append("C:" + nf.format((stat[c][0]+stat[c][1])/all*100) );
			protocol.append("%[S:" + nf.format(stat[c][0]/all*100) );
			protocol.append("%,D:" + nf.format(stat[c][1]/all*100) );
			protocol.append("%],F:" + nf.format(stat[c][2]/all*100) );
			protocol.append("%,M:" + nf.format(stat[c][3]/all*100) );
			protocol.append("%,n:"+a+"\n" );
		}
		protocol.append("\n");
		
		if( regex.length > 0 ) {
			for( int c = 0; c < stat.length; c++ ) {
				if( c <regex.length ) {
					protocol.append("\t"+regex[c] + "\t");
				} else {
					protocol.append("\t"+rem + "\t");
				}
			}
			protocol.append("\n");
		}
		protocol.append( get(stat, all,nf,"Complete BUSCOs (C)",0,1) );
		protocol.append( get(stat, all,nf,"Complete and single-copy BUSCOs (S)",0) );
		protocol.append( get(stat, all,nf, "Complete and duplicated BUSCOs (D)",1) );
		protocol.append( get(stat, all,nf,"Fragmented BUSCOs (F)",2) );
		protocol.append( get(stat, all,nf,"Missing BUSCOs (M)",3) );
		protocol.append( "\nTotal BUSCO groups searched\t"+(int)all+"\n" );
		
		return null;
	}
	
	static void add( int[][] stat, HashSet[] hash ) {
		for( int c = 0; c < hash.length; c++ ) {
			int idx;
			switch( hash[c].size() ) {
				case 0: idx=3; break;
				case 1: idx=0; break;
				default: idx=1; break;
			}
			if( idx <= 1 ) {
				stat[c][idx]++;
				hash[c].clear();
			}			
		}
	}
	
	static String get( int[][] stat, double all, NumberFormat nf, String info, int... index ) {
		StringBuffer sb = new StringBuffer(info);
		for( int c = 0; c < stat.length; c++ ) {
			int sum = 0;
			for( int i = 0; i < index.length; i++ ) {
				sum +=stat[c][index[i]];
			}
			sb.append(get(all, sum, nf));
		}
		sb.append("\n");
		return sb.toString();
	}
	
	static String get( double all, int anz, NumberFormat nf ) {
		return "\t"+ anz + "\t" + nf.format(anz/all*100) + "%";
	}
	
	static int getIndex( String value ) {
		switch( value ) {
		case "Complete": return 0;
		case "Duplicated": return -1;
		case "Fragmented": return 2;
		case "Missing": return 3;
		default: throw new IllegalArgumentException( value );
		}
	}

	@Override
	public ToolParameterSet getToolParameters() {
		try {
			return new ToolParameterSet( getToolName(), 
				new FileParameter("BUSCO", "the BUSCO full table based on transcripts/proteins", "tabular", true, new FileExistsValidator()),
				new FileParameter("IDs", "a table with at leat two columns, the first is the gene ID, the second is the transcript/protein ID. The assignment file from the Extractor can be used or a table can be derived by the user from the gene annotation file (gff,gtf)", "tabular", true, new FileExistsValidator()),
				new ParameterSetContainer("subgenomes", "", new ExpandableParameterSet( 
						new SimpleParameterSet(	
								new SimpleParameter(DataType.STRING,"subgenome","regex for contigs/chromosomes of this subgenome", true )
						),  "subgenomes", "regular expression for subgenome contigs/chromsome names", 0 )
				)
			);
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public String getToolName() {
		return getClass().getSimpleName();
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "recomputes BUSCO statistic for genes";
	}

	@Override
	public String getHelpText() {
		return "This tool can be used to compute BUSCO statistics for genes instead of transcripts."
				+ " Proteins of an annotation file can be extracted with **Exctractor**, Proteins can be used to compute BUSCO statistics with BUSCO."
				+ " The full BUSCO table and the assignment file from the **Extractor** can be used as input for this tool."
				+ " Alternatively, a table can be generated from the annotation file that can be used instead of the assignment file."
				+ MORE;
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public ToolResult[] getTestCases(String path) {
		// TODO Auto-generated method stub
		return null;
	}
}
