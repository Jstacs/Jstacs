package projects.gemoma;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.NumberFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Locale;

import de.jstacs.parameters.FileParameter;
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
	
	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads, String temp) throws Exception {
		FileParameter fp = (FileParameter) parameters.getParameterForName("BUSCO");
		String busco = fp.getValue();
		
		fp = (FileParameter) parameters.getParameterForName("IDs");
		String genTranscript = fp.getValue();
		
		HashMap<String,String> trans2gene = new HashMap<String,String>();
		BufferedReader r = new BufferedReader( new FileReader( genTranscript ) );
		String line;
		while( (line=r.readLine()).charAt(0)=='#' );
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			trans2gene.put(split[1], split[0]);
		}
		r.close();
		
		int[] stat = new int[4];
		double all=0;
		r = new BufferedReader( new FileReader( busco ) );
		String old = null;
		HashSet<String> hash = new HashSet<String>();
		int anz=0;
		while( (line=r.readLine()).charAt(0)=='#' ) {
			if( anz < 3 ) protocol.append(line+"\n");
			anz++;
		}
		protocol.append( "# " + "BUSCORecomputer\n" );
		while( (line=r.readLine()) != null ) {
			String[] split = line.split("\t");
			if( old!= null && !split[0].equals(old) ) {
				stat[hash.size()==1?0:1]++;
				all++;
				hash.clear();
				old=null;
			}
			int v = getIndex(split[1]);
			if( v>= 0 ) {
				stat[v]++;
				all++;
			} else {
				if( old == null ) {
					old=split[0];
				}
				String gene = trans2gene.get(split[2]);
				if( gene == null ) {
					gene = split[2];
					protocol.append("Warning no gene found for transcript: " + gene + "\n");
				}
				hash.add(gene);
			}
		}
		r.close();
		
		NumberFormat nf = NumberFormat.getInstance(Locale.US);
		protocol.append("Complete\t" + nf.format((stat[0]+stat[1])/all) + "\n" );
		protocol.append(" Single\t" + nf.format(stat[0]/all) + "\n" );
		protocol.append(" Duplicated\t" + nf.format(stat[1]/all) + "\n" );
		protocol.append("Fragmented\t" + nf.format(stat[2]/all) + "\n" );
		protocol.append("Missing\t" + nf.format(stat[3]/all) + "\n" );
		
		return null;
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
		return new ToolParameterSet( getToolName(), 
			new FileParameter("BUSCO", "the BUSCO full table based on transcripts/proteins", "tabular", true),
			new FileParameter("IDs", "a table with at leat two columns, the first is the gene ID, the second is the transcript/protein ID. The assignment file from the Extractor can be used or a table can be derived by the user from the gene annotation file (gff,gtf)", "tabular", true)
		);
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
