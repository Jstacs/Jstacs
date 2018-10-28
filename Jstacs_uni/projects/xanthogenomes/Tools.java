package projects.xanthogenomes;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;

import de.jstacs.Singleton;
import de.jstacs.algorithms.alignment.Alignment;
import de.jstacs.algorithms.alignment.Alignment.AlignmentType;
import de.jstacs.algorithms.alignment.PairwiseStringAlignment;
import de.jstacs.algorithms.alignment.cost.AffineCosts;
import de.jstacs.algorithms.alignment.cost.MatrixCosts;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.DataSet;
import de.jstacs.data.EmptyDataSetException;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DNAAlphabetContainer;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.Sequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.ArrayHandler;
import de.jstacs.utils.IntList;
import de.jstacs.utils.PFMComparator;
import de.jstacs.utils.ToolBox;


public class Tools {
	
	public static class Aligner{
		
		public static Aligner DEFAULT;
		
		static{
			try {
				DEFAULT = new Aligner( Aligner.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/BLOSUM62.txt" ), Translator.DEFAULT.getProteinAlphabet() );
			} catch ( Exception e ) { }
		}
		
		private Alignment align;
		private double[][] matrix;
		
		
		public double[][] getMatrix() throws CloneNotSupportedException {
			return ArrayHandler.clone(matrix);
		}

		public Aligner(String blosum, AlphabetContainer prot) throws NumberFormatException, FileNotFoundException, IOException, WrongAlphabetException, CloneNotSupportedException{
			this(new FileInputStream( blosum ),prot);
		}
		
		public Aligner(InputStream blosum, AlphabetContainer prot) throws IOException, NumberFormatException, WrongAlphabetException, CloneNotSupportedException {
			BufferedReader reader = new BufferedReader( new InputStreamReader( blosum ) );
			
			String str = null;
			
			String[] header = null;
			matrix = null;
			while( (str = reader.readLine()) != null ){
				if(!str.startsWith( "#" ) ){
					if(header == null){
						str = str.trim();
						header = str.split("\\s+");
						matrix = new double[header.length][header.length];
					}else{
						String[] pars = str.split("\\s+");
						for(int i=1;i<pars.length;i++){
							//System.out.println(">>"+pars[0].trim()+"<<>>"+header[i-1].trim()+"<<");
							matrix[(int)prot.getCode( 0, pars[0].trim() )][ (int)prot.getCode( 0, header[i-1].trim() ) ] = -Double.parseDouble( pars[i] );
						}
					}
				}
			}
			reader.close();
			
			AffineCosts cost = new AffineCosts(3, new MatrixCosts(matrix, 1));
			align = new Alignment( cost );
		}
		
		public PairwiseStringAlignment align(Sequence tale1, Sequence tale2, AlignmentType type){
			
			PairwiseStringAlignment alignment = align.getAlignment( type, tale1, tale2 );
			
			return alignment;
			
		}
		
	}
	
	public static class ProteinAlphabetContainer extends AlphabetContainer implements Singleton {

		public static final ProteinAlphabetContainer SINGLETON = get();
		
		private static ProteinAlphabetContainer get(){
			ProteinAlphabetContainer res = null;
			try {
				res = new ProteinAlphabetContainer();
			} catch ( Exception doesnothappen ) {
				doesnothappen.printStackTrace();
			}
			
			return res;
		}
		
		/**
		 * @param abc
		 * @throws DoubleSymbolException 
		 * @throws IllegalArgumentException 
		 */
		private ProteinAlphabetContainer( ) throws IllegalArgumentException, DoubleSymbolException {
			super( new DiscreteAlphabet( false, new String[]{"I", "L", "V", "F", "M", "C", "A", "G", "P", "T", "S", "Y", "W", "Q", "N", "H", "E", "D", "K", "R", "*", "B", "Z", "X", "-"} ) );
		}
		
	}
	
	public static class Translator{
		public static Translator DEFAULT;
		static{
			try{
				DEFAULT = new Translator( Translator.class.getClassLoader().getResourceAsStream( "projects/xanthogenomes/genetic_code.txt" )  );
			}catch(Exception e){e.printStackTrace();}
		}
		
		private AlphabetContainer prot;
		private int[][][] codonMap;
		
		public Translator(String genCode) throws IllegalArgumentException, FileNotFoundException, IOException, DoubleSymbolException, WrongAlphabetException{
			this(new FileInputStream( genCode ));
		}
		
		public Translator(InputStream genCode) throws IOException, IllegalArgumentException, DoubleSymbolException, WrongAlphabetException {
			
			AlphabetContainer dna = DNAAlphabetContainer.SINGLETON;
			
			BufferedReader reader = new BufferedReader( new InputStreamReader( genCode ));
			
			String str = null;
			
			HashMap<String,String[]> list = new HashMap<String,String[]>();
			
			while( (str = reader.readLine()) != null ){
				String[] code = str.split( "\t" );
				String[] vals = code[1].split( ",\\s+" );
				list.put( code[0], vals );
			}
			list.put( "B", new String[0] );
			list.put( "Z", new String[0] );
			list.put( "X", new String[0] );
			list.put( "-", new String[0] );
			String[] key = list.keySet().toArray( new String[0] );
			
			prot = ProteinAlphabetContainer.SINGLETON;
			
			codonMap = new int[(int)dna.getAlphabetLengthAt( 0 )][(int)dna.getAlphabetLengthAt( 0 )][(int)dna.getAlphabetLengthAt( 0 )];
			
			for(int i=0;i<key.length;i++){
				int sym = (int) prot.getCode( 0, key[i] );
				String[] vals = list.get( key[i] );
				for(int j=0;j<vals.length;j++){
					Sequence val = Sequence.create(dna,vals[j].trim());
					codonMap[val.discreteVal( 0 )][val.discreteVal(1)][val.discreteVal(2)] = sym;
				}
				
			}
			
			reader.close();
			
		}
		
		public AlphabetContainer getProteinAlphabet(){
			return prot;
		}
		
		public Sequence translate(Sequence dna, int readingFrame) throws WrongAlphabetException, WrongSequenceTypeException{
			int[] prot = new int[(dna.getLength()-readingFrame)/3];
			for(int i=readingFrame;i+2<dna.getLength();i+=3){
				prot[i/3] = codonMap[dna.discreteVal( i )][dna.discreteVal( i+1 )][dna.discreteVal( i+2 )];
			}
			return new IntSequence( this.prot, prot );
		}
		
		
		public static void main(String[] args) throws IllegalArgumentException, WrongAlphabetException, WrongSequenceTypeException, IOException, DoubleSymbolException{
			Sequence dna = Sequence.create( DNAAlphabetContainer.SINGLETON, "ATGGTGCTGATCTAG" );
			
			Translator trans = new Translator( "projects/xcvgenomes/genetic_code.txt" );
			
			System.out.println(trans.translate( dna, 0 ));
			
		}
		
		
	}
	
	public static DataSet extractRVDs(DataSet tales) throws EmptyDataSetException, WrongAlphabetException{
		Sequence[] seqs = new Sequence[tales.getNumberOfElements()];
		for(int i=0;i<tales.getNumberOfElements();i++){
			seqs[i] = tales.getElementAt( i ).getSubSequence( 11, 2 );
		}
		return new DataSet("",seqs);
	}
	
	public static Sequence[] translate(Sequence[] tales, Translator t) throws WrongAlphabetException, WrongSequenceTypeException, EmptyDataSetException, IllegalArgumentException, IOException, DoubleSymbolException{
				
		Sequence [] seqs = new Sequence[tales.length];
		for(int i=0;i<tales.length;i++){
			seqs[i] = t.translate(tales[i],0);
		}
		
		return seqs;
		
	}
	
	
	public static Sequence getConsensusSequence(DataSet alignment) throws WrongAlphabetException, WrongSequenceTypeException{
		
		AlphabetContainer prot = alignment.getAlphabetContainer();
		
		double[][] pfm = PFMComparator.getPFM( alignment );
		
		
		IntList il = new IntList();
		for(int i=0;i<pfm.length;i++){
			int idx = ToolBox.getMaxIndex( pfm[i] );
			if( ! "-".equals(prot.getSymbol( 0, idx ) ) ){
				il.add( idx );
			}
		}
		
		return new IntSequence( prot, il.toArray() );
		
		
	}
	
	
}
