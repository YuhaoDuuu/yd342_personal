package org.borgmoea.spi;

import java.util.Properties;

import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Problem;
import org.moeaframework.core.spi.AlgorithmProvider;
import org.moeaframework.util.TypedProperties;

/**
 * Enables the Borg MOEA to be instantiated within the MOEA Framework.
 */
public class BorgProvider extends AlgorithmProvider {

	@Override
	public Algorithm getAlgorithm(String name, Properties properties, Problem problem) {
		if (name.equalsIgnoreCase("borg")) {
			return new BorgAlgorithm(problem, new TypedProperties(properties));
		}
		
		return null;
	}

}
