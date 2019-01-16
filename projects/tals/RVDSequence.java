package projects.tals;


import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.IntSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.io.SymbolExtractor;

public class RVDSequence extends IntSequence{
	protected AlphabetContainer alphabetsAminoAcidTwelve;
	protected AlphabetContainer alphabetsAminoAcidThirteen;
	private static int[] discreteValRVDTwelve;
	private static int[] discreteValRVDThirteen;
	
	public RVDSequence(AlphabetContainer alphabetsAminoAcidTwelve,AlphabetContainer alphabetsAminoAcidThirteen, String sequence)
			throws WrongAlphabetException, WrongSequenceTypeException, IllegalArgumentException, DoubleSymbolException {
		super(getContainerRVD(alphabetsAminoAcidTwelve,alphabetsAminoAcidThirteen),null, sequence,"-");
		this.alphabetsAminoAcidTwelve=alphabetsAminoAcidTwelve;
		this.alphabetsAminoAcidThirteen=alphabetsAminoAcidThirteen;
		
		
	}
	public int discreteValTwelve( int discreteValRVD ) {

		return discreteValRVDTwelve[discreteValRVD];
	}
	
	public int discreteValThirteen( int discreteValRVD ) {

		return discreteValRVDThirteen[discreteValRVD];
	}
	
	public static AlphabetContainer getContainerRVD(AlphabetContainer alphabetsAminoAcidTwelve, AlphabetContainer alphabetsAminoAcidThirteen) throws IllegalArgumentException, DoubleSymbolException{
		String[] alphabetRVD=new String[(int)alphabetsAminoAcidTwelve.getAlphabetLengthAt(0)*(int)alphabetsAminoAcidThirteen.getAlphabetLengthAt(0)];
		discreteValRVDTwelve=new int[alphabetRVD.length];
		discreteValRVDThirteen=new int[alphabetRVD.length];
		int k=0;
		for(int i=0;i<(int)alphabetsAminoAcidTwelve.getAlphabetLengthAt(0);i++){
			for(int j=0;j<(int)alphabetsAminoAcidThirteen.getAlphabetLengthAt(0);j++){
				alphabetRVD[k]=alphabetsAminoAcidTwelve.getSymbol(0, i)+alphabetsAminoAcidThirteen.getSymbol(0, j);
				discreteValRVDTwelve[k]=i;
				discreteValRVDThirteen[k]=j;
				k++;
			}
		}
	//	System.out.println(alphabetsAminoAcidTwelve.getDelim());
		AlphabetContainer alphabetContainerRVD=new AlphabetContainer(new DiscreteAlphabet(true,alphabetRVD));
		return alphabetContainerRVD;
	}

}
