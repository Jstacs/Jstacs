package projects.dimont;

import de.jstacs.DataType;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;


public class DimontPredictorParameterSet extends SimpleParameterSet {

	public static final String HOME = "home";
	public static final String DATA = "data";
	public static final String INFIX = "infix";
	public static final String CHIPPER = "dimont";
	public static final String PVAL = "p-value";
	public static final String VALUE_TAG = "value";
	public static final String WEIGHTING_FACTOR = "weightingFactor";
	
	public static final String[] PREFIX = {
	                                       CHIPPER, HOME, DATA, INFIX, VALUE_TAG, WEIGHTING_FACTOR, PVAL
	};
	
	public DimontPredictorParameterSet() throws Exception {
		super();
		parameters.add( new SimpleParameter( DataType.STRING, "Dimont", "The file name of the file containing the trained Dimont classifier (absolute path, .xml).", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Home directory", "The path to the directory containing the input file. Output files are written to this directory as well.", true, "./" ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Input file", "The file name of the file containing the input sequences in annotated FastA format (see readme)", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Infix", "a infix to be used for all output files (sequence logos, predicted binding sites)", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Value tag", "The tag for the value information in the FastA-annotation of the input file", true ) );
		parameters.add( new SimpleParameter( DataType.STRING, "Weighting factor", "The value for weighting the data; either a value between 0 and 1, or a description relative to the standard deviation (e.g. +4sd)", true, "" + 0.2 ) );
		
		
		
		parameters.add( new SimpleParameter( DataType.DOUBLE, "p-value", "The maximum p-value allowed for predicted binding sites", true, new NumberValidator<Double>( 0.0, 1.0 ), 1E-3 ) );
		
		
	}
	
}
