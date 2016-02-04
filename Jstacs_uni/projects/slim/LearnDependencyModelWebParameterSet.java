package projects.slim;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;

import de.jstacs.DataType;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.SparseSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.parameters.FileParameter;
import de.jstacs.parameters.ParameterSet;
import de.jstacs.parameters.SelectionParameter;
import de.jstacs.parameters.SimpleParameter;
import de.jstacs.parameters.SimpleParameterSet;
import de.jstacs.parameters.validation.NumberValidator;
import de.jstacs.tools.ui.galaxy.DataColumnParameter;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Pair;


public class LearnDependencyModelWebParameterSet extends ParameterSet {

	public enum ModelType{
		IMM("Inhomogeneous Markov model"),
		BT_EAR("Bayesian tree (EAR)"),
		BT_MI("Bayesian tree (MI)"),
		SLIM("Slim"),
		LSLIM("LSlim");
		
		private String printname;
		
		private ModelType(String printname){
			this.printname = printname;
		}
		
		public String getPrintname(){
			return printname;
		}
		
	}
	
	
	public LearnDependencyModelWebParameterSet() throws Exception {
		
		this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{"Annotated FastA","Tabular"}, new Object[]{
		           new SimpleParameterSet( new FileParameter( "Input sequences", "The input sequences for learning the dependency model (can be uploaded using &quot;GetData&quot; -&gt; &quot;Upload File&quot;), annotated FastA format. The required format is described in the help section.", "fasta", true ) ),
		           new SimpleParameterSet( new FileParameter( "Input sequences", "The input sequences for learning the dependency model (can be uploaded using &quot;GetData&quot; -&gt; &quot;Upload File&quot;), tabular format.", "tabular", true ),
		        		   new DataColumnParameter("Input sequences","Sequence column","The column containing the sequence data",true,1),
		        		   new DataColumnParameter("Input sequences","Signal column","The column containing the signal data for each sequence",true,2)
		        		   )
		}, "Input data", "Select the input data format and set input parameters", true ) );
		
		
		this.parameters.add( new SelectionParameter( DataType.PARAMETERSET, new String[]{ModelType.IMM.getPrintname(),ModelType.BT_EAR.getPrintname(),ModelType.BT_MI.getPrintname(),ModelType.SLIM.getPrintname(),ModelType.LSLIM.getPrintname()}, 
				new Object[]{
				             new SimpleParameterSet( new SimpleParameter(DataType.INT, "Order", "The order of the Markov model", true, new NumberValidator<Integer>(0,2), 1 ) ),
				             new SimpleParameterSet(  ),
				             new SimpleParameterSet(  ),
				             new SimpleParameterSet(  ),
				             new SimpleParameterSet( new SimpleParameter(DataType.INT, "Distance", "The maximum distance of the LSlim model", true, new NumberValidator<Integer>(1,Integer.MAX_VALUE), 5 ) )
				}, "Model type", "Define the type of the dependency model and (if required) additional parameters", true ) );
		
		
		this.parameters.add( new SimpleParameter(DataType.INT, "Background order", "The order of the background model, -1 for uniform distribution", true, new NumberValidator<Integer>(-1,4), -1 ) );
		
		this.parameters.add( new SimpleParameter( DataType.DOUBLE, "Equivalent sample size", "Reflects the strength of the prior on the model parameters.", true, new NumberValidator<Double>( 0.0, Double.MAX_VALUE ), 4.0 ) );
		
		
	}
	
	
	public ModelType getModelType(){
		SelectionParameter sel = (SelectionParameter)this.parameters.get( 1 );
		int selected = sel.getSelected();
		if(selected == 0){
			return ModelType.IMM;
		}else if(selected == 1){
			return ModelType.BT_EAR;
		}else if(selected == 2){
			return ModelType.BT_MI;
		}else if(selected == 3){
			return ModelType.SLIM;
		}else if(selected == 4){
			return ModelType.LSLIM;
		}else{
			throw new RuntimeException( "Unknown model type" );
		}
	}
	
	public double getESS(){
		return (Double)this.parameters.get( 3 ).getValue();
	}
	
	public int getBgOrder(){
		return (Integer)this.parameters.get( 2 ).getValue();
	}
	
	public Integer getOrder(){
		SelectionParameter sel = (SelectionParameter)this.parameters.get( 1 );
		int selected = sel.getSelected();
		if(selected == 0){
			return (Integer)((ParameterSet)sel.getValue()).getParameterAt( 0 ).getValue();
		}else if(selected == 4){
			return -(Integer)((ParameterSet)sel.getValue()).getParameterAt( 0 ).getValue();
		}else{
			return null;
		}
	}
	
	public Pair<DataSet,double[]> getData() throws FileNotFoundException, WrongSequenceTypeException, WrongAlphabetException, EmptyDataSetException, IOException{
		SelectionParameter sel = (SelectionParameter)this.parameters.get( 0 );
		
		if(sel.getSelected() == 0){
			//SplitSequenceAnnotationParser parser = new SplitSequenceAnnotationParser( ":", ";" );
			String filename = (String)((ParameterSet)sel.getValue()).getParameterAt( 0 ).getValue();
			DataSet data= SparseSequence.getDataSet( DNAAlphabetContainer.SINGLETON, new SparseStringExtractor( filename, '>', new SplitSequenceAnnotationParser(":",";") ) );
			
			DoubleList list = new DoubleList();
			
			for(int i=0;i<data.getNumberOfElements();i++){
				list.add( Double.parseDouble( data.getElementAt( i ).getSequenceAnnotationByType( "signal", 0 ).getIdentifier() ) );
			}
			
			return new Pair<DataSet,double[]>(data,list.toArray());
		}else{
			
			ParameterSet parameters = (ParameterSet)sel.getValue();
			
			int seqCol = (Integer)parameters.getParameterAt( 1 ).getValue()-1;
			int valCol = (Integer)parameters.getParameterAt( 2 ).getValue()-1;
			BufferedReader read = new BufferedReader( new FileReader( ( (FileParameter)parameters.getParameterAt( 0 ) ).getFileContents().getFilename() ) );
			LinkedList<Sequence> seqs = new LinkedList<Sequence>();
			DoubleList w = new DoubleList();
			String str = null;
			while( (str = read.readLine()) != null ){
				if(!str.trim().startsWith( "#" )){
					String[] parts = str.split( "\t" );

					seqs.add( Sequence.create( DNAAlphabetContainer.SINGLETON, parts[seqCol], "" ) );
					w.add( Double.parseDouble( parts[valCol] ) );
				}
			}
			
			read.close();
			
			Sequence[] seqs2 = seqs.toArray(new Sequence[0]);
			double[] w2 = w.toArray();
			
			return new Pair<DataSet, double[]>( new DataSet("",seqs2), w2 );
			
		}
		
	}
	
	
	
}
