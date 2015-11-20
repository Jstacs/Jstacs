package projects.dimont;

import de.jstacs.DataType;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;

/**
 * This class is a container for all parameters of Dimont. It also parses the parameter from Strings.
 *  
 * @author Jens Keilwagen
 */
public class DimontParameterSet extends SimpleParameterSet {
	
	public static final String HOME = "home";
	public static final String DATA = "data";
	public static final String INFIX = "infix";
	public static final String LENGTH = "motifWidth";
	public static final String STARTS = "starts";
	public static final String MOTIF_ORDER = "motifOrder";
	public static final String BG_ORDER = "bgOrder";
	public static final String POSITION_TAG = "position";
	public static final String VALUE_TAG = "value";
	public static final String SD = "sd";
	public static final String WEIGHTING_FACTOR = "weightingFactor";
	public static final String ESS = "ess";
	public static final String DELETE = "delete";
	public static final String THREADS = "threads";

	public static final String[] PREFIX = {
        HOME, DATA, INFIX, POSITION_TAG, VALUE_TAG, SD, WEIGHTING_FACTOR, STARTS, LENGTH, MOTIF_ORDER, BG_ORDER, ESS, DELETE, THREADS
    };
	
	public DimontParameterSet() throws Exception {
		super();

		parameters.add( new SimpleParameter( DataType.STRING, "Home directory", "The path to the directory containing the input file. Output files are written to this directory as well.", true, "./" ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Input file", "The file name of the file containing the input sequences in annotated FastA format (see readme)", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Infix", "a infix to be used for all output files (model, sequence logos, predicted binding sites)", true ) );
		
		parameters.add( new SimpleParameter( DataType.STRING, "Position tag", "The tag for the position information in the FastA-annotation of the input file", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Value tag", "The tag for the value information in the FastA-annotation of the input file", true ) );
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Standard deviation", "The standard deviation of the position distribution centered at the position specified by the position tag", true, new NumberValidator<Double>( 1.0, 1E4 ), 75.0 ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Weighting factor", "The value for weighting the data; either a value between 0 and 1, or a description relative to the standard deviation (e.g. +4sd)", true, "" + 0.2 ) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Starts", "The number of pre-optimization runs.", true, new NumberValidator<Integer>(1,100), 20 ) );
		
		
		parameters.add( new SimpleParameter( DataType.INT, "Initial motif width", "The motif width that is used initially, may be adjusted during optimization.", true, new NumberValidator<Integer>(1,50), 15 ) );
		
		parameters.add( new SimpleParameter( DataType.INT, "Markov order of motif model", "The Markov order of the model for the motif.", true, new NumberValidator<Integer>(0,3), 0 ) );
		parameters.add( new SimpleParameter( DataType.INT, "Markov order of background model", "The Markov order of the model for the background sequence and the background sequence, -1 defines uniform distribution.", true, new NumberValidator<Integer>(-1,5), -1 ) );
		
		
		
		parameters.add( new SimpleParameter( DataType.DOUBLE, "Equivalent sample size", "Reflects the strength of the prior on the model parameters.", true, new NumberValidator<Double>(0d, Double.POSITIVE_INFINITY), 4d ) );
		
		parameters.add( new SimpleParameter( DataType.BOOLEAN, "Delete BSs from profile", "A switch for deleting binding site positions of discovered motifs from the profile before searching for futher motifs.", true, true ) );
		
		parameters.add( new SimpleParameter( DataType.INT, "Compute threads", "The number of threads that are use to evaluate the objective function and its gradient.", false, new NumberValidator<Integer>(1,128) ) );
	}
}
