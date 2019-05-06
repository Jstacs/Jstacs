package projects.tals.training;

import java.util.HashMap;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.WrongAlphabetException;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.ReferenceSequenceAnnotation;

public class CachingRVDReferenceSequenceAnnotationParser extends RVDReferenceSequenceAnnotationParser {

	private HashMap<String, ReferenceSequenceAnnotation> refCache;
	
	public CachingRVDReferenceSequenceAnnotationParser(String key, String keyValueDelimiter, String annotationDelimiter,
			AlphabetContainer alphabet12, AlphabetContainer alphabet13) throws IllegalArgumentException {
		super(key, keyValueDelimiter, annotationDelimiter, alphabet12, alphabet13);
		this.refCache = new HashMap<String, ReferenceSequenceAnnotation>();
	}

	@Override
	protected ReferenceSequenceAnnotation getSequenceAnnotation(String seqString)
			throws IllegalArgumentException, WrongAlphabetException, WrongSequenceTypeException, DoubleSymbolException {
		if(refCache.containsKey(seqString)){
			return refCache.get(seqString);
		}else{
			ReferenceSequenceAnnotation annotation = super.getSequenceAnnotation(seqString);
			refCache.put(seqString, annotation);
			return annotation;
		}
	}
	
	

}
