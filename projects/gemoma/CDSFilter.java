package projects.gemoma;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import de.jstacs.tools.Protocol;
import de.jstacs.tools.ui.cli.CLI.SysProtocol;

/**
 * Filters a GFF for CDS and modifies transcripts and genes borders accordingly.
 * 
 * @author Jens Keilwagen
 */
public class CDSFilter {

	public static void main(String[] args) throws IOException {
		reformat(args[0], new SysProtocol() );
	}

	private static String par = "Parent=";
	
	public static void reformat( String input, Protocol protocol ) throws IOException {
		BufferedReader r;
		String line;
		String[] split=null;
		int idx, h;
		
		//read transcripts
		r = Tools.openGzOrPlain(input);
		HashMap<String, HashSet<String>> geneTrans = new HashMap<String, HashSet<String>>();
		HashMap<String, ArrayList<String[]>> transCDS = new HashMap<String, ArrayList<String[]>>();
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
			
			switch( split[2] ) {
				case "CDS":
					if( split[8].indexOf(par)<0 ) {
						idx = split[8].indexOf("ID=");
						if( idx>= 0 ) {
							split[8] = split[8].replace("ID=", par);
						} else {
							throw new IllegalArgumentException("The line does not contain ID or Parent: " + line);
						}
					}
					idx = split[8].indexOf(par) + par.length();
					h = split[8].indexOf(';',idx);
					String[] tra = split[8].substring(idx, h).split(",");
					for( String tr: tra ) {
						addCDSPart( tr, split, transCDS );
					}
					break;
				case "mRNA": case "transcript": case "prediction": //XXX
					idx = split[8].indexOf("ID=")+3;
					h = split[8].indexOf(';',idx);
					String tr = split[8].substring(idx, h>0?h:split[8].length() );
					idx = split[8].indexOf(par);
					if( idx>=0 ) {
						idx+=par.length();
						h = split[8].indexOf(';',idx);
					}
					String geneID = idx<0 ? tr+".gene" : split[8].substring(idx, h>0?h:split[8].length() );
					addTranscript( geneID, tr, geneTrans );
					break;
				default:
					//skip
				}
		}
		r.close();

		//write
		Iterator<String> genes = geneTrans.keySet().iterator();
		StringBuffer gOut = new StringBuffer();
		while( genes.hasNext() ) {
			String geneId = genes.next();
			gOut.delete(0, gOut.length());
			int gStart=Integer.MAX_VALUE, gEnd = 0;
			Iterator<String> trans = geneTrans.get(geneId).iterator();
			StringBuffer tOut = new StringBuffer();
			while( trans.hasNext() ) {
				String transId = trans.next();
				tOut.delete(0, tOut.length());
				int start=Integer.MAX_VALUE, end = 0;
				ArrayList<String[]> cds = transCDS.get(transId);
				if( cds.size()>0 ) {
					for( int i = 0; i < cds.size(); i++ ) {
						split = cds.get(i);
						for( int j = 0; j < split.length; j++ ) {
							tOut.append( split[j] + (j+1==split.length?"\n":"\t") );
						}
						start = Math.min(start, Integer.parseInt(split[3]) );
						end = Math.max(end, Integer.parseInt(split[4]) );
					}
					gOut.append( split[0] + "\tderived\ttranscript\t" + start+ "\t" + end +"\t.\t" + split[6] + "\t.\tID="+transId+";Parent="+geneId+"\n");
					gOut.append( tOut );
					gStart = Math.min(start, gStart );
					gEnd = Math.max(end, gEnd );
				}
			}
			System.out.print( split[0] + "\tderived\tgene\t" + gStart+ "\t" + gEnd +"\t.\t" + split[6] + "\t.\tID="+geneId+"\n");			
			System.out.print( gOut );
		}
	}

	public static void addCDSPart( String transID, String[] split, HashMap<String,ArrayList<String[]>> transCDS ) {
		ArrayList<String[]> cds = transCDS.get(transID);
		if( cds == null ) {
			cds = new ArrayList<String[]>();
			transCDS.put(transID, cds);
		}
		cds.add(split);
	}
	
	public static void addTranscript( String geneID, String transID, HashMap<String,HashSet<String>> geneTrans ) {
		HashSet<String> trans = geneTrans.get(geneID);
		if( trans == null ) {
			trans = new HashSet<String>();
			geneTrans.put( geneID, trans );
		}
		trans.add(transID);
	}
}
