package projects.gemoma;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import de.jstacs.DataType;
import de.jstacs.parameters.ExpandableParameterSet;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.Parameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.ParameterSetContainer;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.tools.JstacsTool;
import de.jstacs.tools.ProgressUpdater;
import de.jstacs.tools.Protocol;
import de.jstacs.tools.ToolParameterSet;
import de.jstacs.tools.ToolResult;

/**
 * The class creates statistics for introns.
 *  
 * @author Jens Keilwagen
 */
public class CheckIntrons implements JstacsTool {

	@Override
	public ToolParameterSet getToolParameters() {
		try{
			return new ToolParameterSet( getShortName(), 
					new FileParameter( "target genome", "The target genome file (FASTA), i.e., the target sequences in the blast run. Should be in IUPAC code", "fasta", true ),					
					new ParameterSetContainer( "introns", "", new ExpandableParameterSet( new SimpleParameterSet(	
							new FileParameter( "introns", "Introns (GFF), which might be obtained from RNA-seq", "gff", false )
						), "introns", "", 1 ) ),
					
					new SimpleParameter( DataType.BOOLEAN, "verbose", "A flag which allows to output a wealth of additional information per transcript", true, false )
			);
		}catch(Exception e){
			e.printStackTrace();
			throw new RuntimeException();
		}
	}

	@Override
	public ToolResult run(ToolParameterSet parameters, Protocol protocol, ProgressUpdater progress, int threads) throws Exception {
		boolean verbose = (Boolean) parameters.getParameterForName("verbose").getValue();
		String targetGenome = (String) parameters.getParameterForName("target genome").getValue(); 
		int reads = 1;
		ExpandableParameterSet introns = (ExpandableParameterSet)((ParameterSetContainer)parameters.getParameterAt(1)).getValue();
		
		HashMap<String,String> seqs = Tools.getFasta(targetGenome,20,' ');
		
		ArrayList<String> fName = new ArrayList<String>();
		for( int i = 0; i < introns.getNumberOfParameters(); i++ ) {
			Parameter y = ((ParameterSet)introns.getParameterAt(i).getValue()).getParameterAt(0);
			if( y.isSet() ) {
				fName.add(y.getValue().toString());
			}
		}
		if( fName.size()>0 ) {
			HashMap<String,int[]> diNucl = new HashMap<String,int[]>();
			HashMap<String, int[][][]>[] res = GeMoMa.readIntrons( reads, protocol, verbose, seqs, diNucl, fName.toArray(new String[fName.size()]) );
			String[] keys = diNucl.keySet().toArray(new String[0]);
			Arrays.sort( keys );
			for( int i = 0; i < keys.length; i++ ) {
				int[] v = diNucl.get( keys[i] );
				protocol.append( keys[i] + "\t" + v[0] + "\t" + v[1] + "\n" );
			}
		}
		return null;
	}

	@Override
	public String getToolName() {
		return "CheckIntrons";
	}

	@Override
	public String getToolVersion() {
		return GeMoMa.VERSION;
	}

	@Override
	public String getShortName() {
		return getToolName();
	}

	@Override
	public String getDescription() {
		return "creates statistics for introns";
	}

	@Override
	public String getHelpText() {
		return "The tool checks the distribution of introns on the strands and the dinucleotide distribution at splice sites.";
	}

	@Override
	public ResultEntry[] getDefaultResultInfos() {
		return null;
	}
}
