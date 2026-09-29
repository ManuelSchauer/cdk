/*
 * CDK-Descriptor-Calculation
 * Copyright (C) 2026 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * Source code is available at <https://github.com/JonasSchaub/CDK-Descriptor-Calculation>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.openscience.cdk.qsar;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.qsar.descriptors.molecular.ALOGPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.APolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AcidicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AminoAcidCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticAtomsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AromaticBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability;
import org.openscience.cdk.qsar.descriptors.molecular.BCUTDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BPolDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BasicGroupCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.BondCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CarbonTypesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathClusterDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ChiPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.EccentricConnectivityIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FMFDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalCSP3Descriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FractionalPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.FragmentComplexityDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.HybridizationRatioDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.JPlogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KappaShapeIndicesDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.KierHallSmartsDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LargestPiSystemDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LongestAliphaticChainDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MDEDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MannholdLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanNumberDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RotatableBondsCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.RuleOfFiveDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SmallRingDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.SpiroAtomCountDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.TPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VABCDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.VAdjMaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.HydrogenState;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

/**
 * Descriptor related calculations based on the CDK for the enrichment of data vectors.
 * <p>
 *     The available descriptors and fingerprints, including their characteristics (e.g. the number of
 *     calculated components), are defined by the {@link Descriptor} enum.
 * </p>
 * <p>
 *     There are 4 different "public static boolean setDescriptorsForMolecules...()" methods with different forms of
 *     (parallelized) calculation. There is single molecule processing or batch processing for large datastructures.
 *     Single and Batch processing can be used with either an IAtomContainer array or with a String Array of SMILES codes.
 * </p>
 * <p>
 *     The four primary public calculation methods differ along two independent axes – the <i>input type</i>
 *     (pre-parsed {@link IAtomContainer} objects vs. raw SMILES strings) and the
 *     <i>parallelization strategy</i> (molecule-level parallelization vs. batch-level parallelization) – and are
 *     summarized in the table below. Every method accepts a {@code boolean isParallelCalculation} flag; passing
 *     {@code false} disables parallelization and the method runs sequentially instead.
 *     <br><br>
 *     <table border="1">
 *         <caption>Overview of the four public batch-calculation methods</caption>
 *         <tr>
 *             <th>Method</th>
 *             <th>Input type</th>
 *             <th>Parallelisation unit</th>
 *             <th>Extra parameter</th>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesByMoleculeParallelization}</td>
 *             <td>{@link IAtomContainer IAtomContainer[]}</td>
 *             <td>One thread per <b>molecule</b> (Java parallel stream over all molecule indices)</td>
 *             <td>–</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesByBatchParallelization}</td>
 *             <td>{@link IAtomContainer IAtomContainer[]}</td>
 *             <td>One thread per <b>batch</b> of molecules (parallel stream over batch indices)</td>
 *             <td>{@code batchSize} – number of molecules per batch</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculesBySmilesStringParallelization}</td>
 *             <td>{@code String[]} (SMILES)</td>
 *             <td>One thread per <b>molecule</b> (Java parallel stream over all molecule indices)</td>
 *             <td>{@code electronDonationModel} – aromaticity model applied during SMILES parsing</td>
 *         </tr>
 *         <tr>
 *             <td>{@link #setDescriptorsForMoleculeBySmilesStringsBatchParallelization}</td>
 *             <td>{@code String[]} (SMILES)</td>
 *             <td>One thread per <b>batch</b> of molecules (parallel stream over batch indices)</td>
 *             <td>{@code batchSize}, {@code electronDonationModel}</td>
 *         </tr>
 *     </table>
 *     <br>
 *     <b>Molecule-level parallelization</b> ({@code ...ByMoleculeParallelization}) dispatches every molecule as an
 *     independent work item to Java's common fork-join pool via {@link IntStream#parallel()}.
 *     This provides fine-grained load balancing and is generally the preferred choice when the molecule set is
 *     large and individual descriptor calculations vary in cost.
 *     <br><br>
 *     <b>Batch-level parallelization</b> ({@code ...ByBatchParallelization}) first divides the molecule array into
 *     fixed-size batches of {@code batchSize} molecules and then submits one batch per fork-join task. Within
 *     each batch, molecules are processed sequentially. This strategy can reduce fork-join overhead for very fast
 *     descriptors or when the overhead of spawning one task per molecule would dominate.
 *     <br><br>
 *     <b>SMILES-based methods</b> additionally parse each SMILES string into an
 *     {@link IAtomContainer} and apply aromaticity perception (using the supplied
 *     {@link ElectronDonation} model, or {@code Daylight} by default if
 *     {@code null} is passed) before descriptor calculation. Pre-parsed {@code IAtomContainer} methods skip this
 *     step; aromaticity must have been perceived on the containers beforehand (e.g. via
 *     {@link #setAromaticity(IAtomContainer, ElectronDonation)}).
 *     <br><br>
 *     All four methods share the same result convention: they return {@code true} if every descriptor value was
 *     computed without producing {@link Float#NaN}, and {@code false} if at least one NaN was encountered. NaN
 *     positions are additionally recorded in the caller-supplied {@code nanPositionsList} list as
 *     {@code int[]{moleculeIndex, componentIndex}} pairs. For parallel executions this list <b>must</b> be
 *     thread-safe (e.g. {@code Collections.synchronizedList(new LinkedList<>())}).
 *     <br><br>
 *     Thread safety for fingerprints is handled internally via a pool of
 *     {@link IFingerprinter} instances (one pool per fingerprint descriptor),
 *     because CDK fingerprinter objects are not thread-safe. The pool size defaults to 4 and can be adjusted
 *     with {@link #setFingerprintPoolSize(int)}.
 * <p>
 *     Example usage of the DescriptorCalculator class:
 * </p>
 * <pre>{@code
 * // Preprocessing
 * // A: Molecule preprocessing:
 * // 1: Create a molecule
 * IAtomContainer molecule = ...;
 *
 * // 2. IMPORTANT: Perform aromaticity perception before descriptor calculation
 * DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
 *
 * // 3. Create a molecule array
 * IAtomContainer[] molecules = new IAtomContainer[] { molecule };
 *
 * // B: String preprocessing
 * // 1: Create a String array of SMILES codes
 * String[] moleculeSmilesStrings = new String[] {"CCO", "CCC(O)O", ...}
 *
 * // Define parameters
 * // 1: Define descriptor array
 * Descriptor[] descriptors = new Descriptor[] {
 *     Descriptor.MOLECULAR_WEIGHT,
 *     Descriptor.TPSA,
 *     Descriptor.A_LOG_P
 * };
 * // or use a method to get specific descriptors:
 * Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(false, false, false, false); // e.g., user-defined selection
 *
 * // 2. Determine total number of components needed (because some descriptors return multiple values, it is not equal to the number of descriptors)
 * int componentCount = Descriptor.getNumberOfComponents(descriptors);
 *
 * // 3. Determine total number of molecules
 * int moleculeCount = molecules.length;
 *
 * // 4. Create data matrix for results with appropriate size (molecules as rows, descriptors as columns)
 * float[][] matrix = new float[moleculeCount][componentCount];
 *
 * // 5. Set startIndex to 0 for filling the matrix with descriptor values starting at the 0th column
 * int startIndex = 0; // if the matrix is part of a larger data structure, set startIndex accordingly
 *
 * // 6. Choose whether to use parallel calculation
 * boolean isParallelCalculation = true; // false for sequential calculation
 *
 * // 7. Create a synchronized List for keeping track of NaN positions
 * List<int[]> nanPositionsListSynchronized = Collections.synchronizedList(new LinkedList<>());
 *
 * // For Batch Processing:
 * // Define batch size
 * int batchSize = 100;
 *
 * // For SMILES Processing
 * // Define anElectronDonation model for aromaticity handling
 * ElectronDonation electronDonationModel = Aromaticity.Model.Daylight;
 *
 * // Calculate descriptors (choose one of the following methods)
 *
 * // Option A: High-Performance-Parallel Descriptor Calculation with IAtomContainer array
 * boolean result = DescriptorCalculator.setDescriptorsForMoleculesByMoleculeParallelization(
 *                      descriptors,
 *                      molecules,
 *                      matrix,
 *                      startIndex,
 *                      isParallelCalculation,
 *                      nanPositionsListSynchronized
 * );
 * // you can also use DescriptorCalculator.setDescriptorsForMoleculesByBatchParallelization(...) -> batchSize needs to be specified accordingly
 *
 * // Option B: High-Performance-Parallel Descriptor Calculation with SMILES string array
 * boolean result = DescriptorCalculator.setDescriptorsForMoleculesBySmilesStringParallelization(
 *                      descriptors,
 *                      moleculeSmilesStrings,
 *                      matrix,
 *                      startIndex,
 *                      electronDonationModel,
 *                      isParallelCalculation,
 *                      nanPositionsListSynchronized
 * );
 * // you can also use DescriptorCalculator.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(...) -> batchSize needs to be specified accordingly
 * }</pre>
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 * @version 1.0.0
 */
public final class DescriptorCalculator {
    /**
     * Logger of this class.
     */
    static final ILoggingTool LOGGER = LoggingToolFactory.createLoggingTool(DescriptorCalculator.class);

    /**
     * SMILES Parser for SMILES String batch processing.
     */
    private static final SmilesParser SMILES_PARSER = new SmilesParser(SilentChemObjectBuilder.getInstance());

    /**
     * EnumMap that maps a descriptor enum constant to an instance of its associated CDK descriptor class.
     */
    private static final EnumMap<Descriptor, IMolecularDescriptor> descriptorToCdkObjectMap = new EnumMap<>(Descriptor.class);
    /**
     * EnumMap that maps a fingerprint descriptor to its pool of IFingerprinter instances.
     * Pool size is configurable via {@link #setFingerprintPoolSize(int)} method.
     * Uses BlockingQueue to ensure thread-safe access to fingerprinter instances.
     */
    private static final EnumMap<Descriptor, BlockingQueue<IFingerprinter>> fingerprintPoolMap = new EnumMap<>(Descriptor.class);
    /**
     * Pool size for fingerprinter instances. Default value is 4 which should be sufficient for regular users.
     * Can be changed via the synchronized {@link #setFingerprintPoolSize(int)} method. Note that the variable
     * itself is not thread-safe!
     */
    private static volatile int fingerprintPoolSize = 4;
    /**
     * Size used for all circular fingerprint "descriptors". Default value is 1024.
     * Can be changed via the synchronized {@link #setCircularFingerprintSize(int)} method.
     */
    private static volatile int circularFingerprintSize = Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE;
    /*
     * Static initializer block to populate the descriptorToCdkObjectMap and initialize the fingerprint pool.
     * Note: We use the map and initialize it here (instead of giving each descriptor constant an instance field)
     * to be able to do error handling, and because IFingerprinter and IMolecularDescriptor are different object types.
     */
    static {
        try {
            // MOLECULAR_WEIGHT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.MOLECULAR_WEIGHT, new WeightDescriptor());

            // WIENER_NUMBER
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.WIENER_NUMBER, new WienerNumbersDescriptor());

            // ATOM_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT, new AtomCountDescriptor());

            // ATOM_COUNT_HEAVY
            AtomCountDescriptor atomCountHeavyDescriptor = new AtomCountDescriptor();
            atomCountHeavyDescriptor.setParameters(new Object[] {"#"}); // set heavy atom count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_HEAVY, atomCountHeavyDescriptor);

            // ATOM_COUNT_C
            AtomCountDescriptor atomCountCDescriptor = new AtomCountDescriptor();
            atomCountCDescriptor.setParameters(new Object[] {"C"}); // set carbon as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_C, atomCountCDescriptor);

            // ATOM_COUNT_H
            AtomCountDescriptor atomCountHDescriptor = new AtomCountDescriptor();
            atomCountHDescriptor.setParameters(new Object[] {"H"}); // set hydrogen as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_H, atomCountHDescriptor);

            // ATOM_COUNT_N
            AtomCountDescriptor atomCountNDescriptor = new AtomCountDescriptor();
            atomCountNDescriptor.setParameters(new Object[] {"N"}); // set nitrogen as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_N, atomCountNDescriptor);

            // ATOM_COUNT_O
            AtomCountDescriptor atomCountODescriptor = new AtomCountDescriptor();
            atomCountODescriptor.setParameters(new Object[] {"O"}); // set oxygen as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_O, atomCountODescriptor);

            // ATOM_COUNT_S
            AtomCountDescriptor atomCountSDescriptor = new AtomCountDescriptor();
            atomCountSDescriptor.setParameters(new Object[] {"S"}); // set sulfur as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_S, atomCountSDescriptor);

            // ATOM_COUNT_P
            AtomCountDescriptor atomCountPDescriptor = new AtomCountDescriptor();
            atomCountPDescriptor.setParameters(new Object[] {"P"}); // set phosphorus as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_P, atomCountPDescriptor);

            // ATOM_COUNT_F
            AtomCountDescriptor atomCountFDescriptor = new AtomCountDescriptor();
            atomCountFDescriptor.setParameters(new Object[] {"F"}); // set fluorine as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_F, atomCountFDescriptor);

            // ATOM_COUNT_BR
            AtomCountDescriptor atomCountBrDescriptor = new AtomCountDescriptor();
            atomCountBrDescriptor.setParameters(new Object[] {"Br"}); // set bromine as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_BR, atomCountBrDescriptor);

            // ATOM_COUNT_CL
            AtomCountDescriptor atomCountClDescriptor = new AtomCountDescriptor();
            atomCountClDescriptor.setParameters(new Object[] {"Cl"}); // set chlorine as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_CL, atomCountClDescriptor);

            // ATOM_COUNT_I
            AtomCountDescriptor atomCountIDescriptor = new AtomCountDescriptor();
            atomCountIDescriptor.setParameters(new Object[] {"I"}); // set iodine as the atom to count
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_I, atomCountIDescriptor);

            // H_BOND_ACCEPTOR_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.H_BOND_ACCEPTOR_COUNT, new HBondAcceptorCountDescriptor());

            // H_BOND_DONOR_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.H_BOND_DONOR_COUNT, new HBondDonorCountDescriptor());

            // TPSA
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.TPSA, new TPSADescriptor());

            // LARGEST_CHAIN
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.LARGEST_CHAIN, new LargestChainDescriptor());

            // LONGEST_ALIPHATIC_CHAIN
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.LONGEST_ALIPHATIC_CHAIN, new LongestAliphaticChainDescriptor());

            // MANNHOLD_LOGP
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.MANNHOLD_LOGP, new MannholdLogPDescriptor());

            // BCUT
            BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
            bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BCUT, bcutDescriptor);

            // BOND_COUNT_ALL
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_ALL, new BondCountDescriptor());

            // BOND_COUNT_SINGLE
            BondCountDescriptor bondCountSingleDescriptor = new BondCountDescriptor();
            bondCountSingleDescriptor.setParameters(new Object[]{"s"});
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_SINGLE, bondCountSingleDescriptor);

            // BOND_COUNT_DOUBLE
            BondCountDescriptor bondCountDoubleDescriptor = new BondCountDescriptor();
            bondCountDoubleDescriptor.setParameters(new Object[]{"d"});
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_DOUBLE, bondCountDoubleDescriptor);

            // BOND_COUNT_TRIPLE
            BondCountDescriptor bondCountTripleDescriptor = new BondCountDescriptor();
            bondCountTripleDescriptor.setParameters(new Object[]{"t"});
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_TRIPLE, bondCountTripleDescriptor);

            // B_POL
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.B_POL, new BPolDescriptor());

            // RULE_OF_FIVE
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.RULE_OF_FIVE, new RuleOfFiveDescriptor());

            // AROMATIC_ATOMS_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AROMATIC_ATOMS_COUNT, new AromaticAtomsCountDescriptor());

            // AROMATIC_BONDS_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AROMATIC_BONDS_COUNT, new AromaticBondsCountDescriptor());

            // ROTATABLE_BONDS_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ROTATABLE_BONDS_COUNT, new RotatableBondsCountDescriptor());

            // FMF
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.FMF, new FMFDescriptor());

            // FRACTIONAL_CSP3
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.FRACTIONAL_CSP3, new FractionalCSP3Descriptor());

            // HYBRIDIZATION_RATIO
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.HYBRIDIZATION_RATIO, new HybridizationRatioDescriptor());

            // KAPPA_SHAPE_INDICES
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.KAPPA_SHAPE_INDICES, new KappaShapeIndicesDescriptor());

            // PETITJEAN_NUMBER
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.PETITJEAN_NUMBER, new PetitjeanNumberDescriptor());

            // SPIRO_ATOM_COUNT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.SPIRO_ATOM_COUNT, new SpiroAtomCountDescriptor());

            // V_ADJ_MAT
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.V_ADJ_MAT, new VAdjMaDescriptor());

            // WEIGHTED_PATH
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.WEIGHTED_PATH, new WeightedPathDescriptor());

            // ZAGREB_INDEX
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ZAGREB_INDEX, new ZagrebIndexDescriptor());

            // CARBON_TYPES
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CARBON_TYPES, new CarbonTypesDescriptor());

            // A_LOG_P
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.A_LOG_P, new ALOGPDescriptor());

            // X_LOG_P
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.X_LOG_P, new XLogPDescriptor());

            // JP_LOG_P
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.JP_LOG_P, new JPlogPDescriptor());

            // A_POL
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.A_POL, new APolDescriptor());

            // AUTOCORRELATION_CHARGE
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_CHARGE, new AutocorrelationDescriptorCharge());

            // AUTOCORRELATION_MASS
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_MASS, new AutocorrelationDescriptorMass());

            // AUTOCORRELATION_POLARIZABILITY
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_POLARIZABILITY, new AutocorrelationDescriptorPolarizability());

            // FRAGMENT_COMPLEXITY
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.FRAGMENT_COMPLEXITY, new FragmentComplexityDescriptor());

            // CHI_CHAIN
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CHI_CHAIN, new ChiChainDescriptor());

            // CHI_CLUSTER
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CHI_CLUSTER, new ChiClusterDescriptor());

            // CHI_PATH_CLUSTER
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CHI_PATH_CLUSTER, new ChiPathClusterDescriptor());

            // CHI_PATH
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CHI_PATH, new ChiPathDescriptor());

            // FRACTIONAL_PSA
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.FRACTIONAL_PSA, new FractionalPSADescriptor());

            // LARGEST_PI_SYSTEM
            LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
            largestPiSystemDescriptor.setParameters(new Object[] {false}); //do not check aromaticity again
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.LARGEST_PI_SYSTEM, largestPiSystemDescriptor);

            // SMALL_RING
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.SMALL_RING, new SmallRingDescriptor());

            // Basic Group Count
            BasicGroupCountDescriptor basicGroupCountDescriptor = new BasicGroupCountDescriptor();
            basicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.BASIC_GROUP_COUNT, basicGroupCountDescriptor);

            // Acidic Group Count
            AcidicGroupCountDescriptor acidicGroupCountDescriptor = new AcidicGroupCountDescriptor();
            acidicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ACIDIC_GROUP_COUNT, acidicGroupCountDescriptor);

            // AminoAcidCount
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.AMINO_ACID_COUNT, new AminoAcidCountDescriptor());

            // Kier-Hall SMARTS
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.KIER_HALL_SMARTS, new KierHallSmartsDescriptor());

            // ECCENTRIC_CONNECTIVITY_INDEX
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.ECCENTRIC_CONNECTIVITY_INDEX, new EccentricConnectivityIndexDescriptor());

            // MDE
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.MDE, new MDEDescriptor());

            // VABC
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.VABC, new VABCDescriptor());

            // CPSA
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.CPSA, new CPSADescriptor());

            // GRAVITATIONAL_INDEX
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.GRAVITATIONAL_INDEX, new GravitationalIndexDescriptor());

            // MOMENT_OF_INERTIA
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.MOMENT_OF_INERTIA, new MomentOfInertiaDescriptor());

            // WHIM
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.WHIM, new WHIMDescriptor());

            // LENGTH_OVER_BREADTH
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.LENGTH_OVER_BREADTH, new LengthOverBreadthDescriptor());

            // PETITJEAN_SHAPE_INDEX
            DescriptorCalculator.descriptorToCdkObjectMap.put(Descriptor.PETITJEAN_SHAPE_INDEX, new PetitjeanShapeIndexDescriptor());

            // Add new descriptor information here!
            // (The enum constant is added in Descriptor, same marker. For a fingerprint, add a pool in
            // initializeFingerprintPools() instead of an entry here.)

            // Initialize fingerprint pools with configurable pool size (default: 4)
            DescriptorCalculator.initializeFingerprintPools();

        } catch (Exception exception) {
            throw new UnsupportedOperationException("Failed to initialize descriptors, this should never happen. ", exception);
        }
    }

    /**
     * Private constructor, this class only provides static methods and is not meant to be instantiated.
     */
    private DescriptorCalculator() {
        // not instantiable
    }

    /**
     * Calculates descriptor or fingerprint values for a molecule and stores them in the result vector.
     * <p>
     * This is the core calculation method that dispatches to either fingerprint or molecular descriptor
     * calculation based on the descriptor type. It retrieves pre-initialized CDK descriptor instances
     * from a shared pool to avoid repeated instantiation overhead.
     * <p>
     * <b>Important:</b> This method does NOT perform input validation. All necessary checks
     * must be performed by the calling public methods before invoking this method.
     * <p>
     * <b>Thread Safety:</b> For fingerprints, this method uses a blocking queue pool
     * ({@link #calculateFingerprintFromPool}) because fingerprinter instances are not thread-safe.
     * For molecular descriptors, shared instances from {@link #descriptorToCdkObjectMap} are used.
     * <p>
     * <b>Result Handling:</b> The method handles four CDK result types:
     * <ul>
     *   <li>{@link DoubleResult} - Single double value (e.g., molecular weight)</li>
     *   <li>{@link IntegerResult} - Single integer value (e.g., atom count)</li>
     *   <li>{@link DoubleArrayResult} - Array of doubles (e.g., BCUT returns 6 values)</li>
     *   <li>{@link IntegerArrayResult} - Array of integers (e.g., amino acid counts)</li>
     * </ul>
     *
     * @param descriptor    the descriptor to calculate (IS NOT CHANGED)
     * @param atomContainer the molecule to calculate descriptors for (IS NOT CHANGED);
     *                        must have aromaticity already perceived if required by the descriptor
     * @param vector         the result vector to store calculated values (MAY BE CHANGED);
     *                        values are stored starting at startIndex
     * @param startIndex     the starting index in vector where results should be written;
     *                        subsequent values are written to startIndex + 1, startIndex + 2, etc.
     * @return Returns true if the calculation was successful and results were stored in the vector, false otherwise.
     */
    private static boolean calculate(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex) throws InterruptedException {
        boolean success = false;
        try {
            if (descriptor.isFingerprint()) {
                return DescriptorCalculator.calculateFingerprintFromPool(descriptor, atomContainer, vector, startIndex);
            } else {
                IMolecularDescriptor cdkDescriptor = DescriptorCalculator.descriptorToCdkObjectMap.get(descriptor);
                if (cdkDescriptor != null) {

                    IDescriptorResult result = cdkDescriptor.calculate(atomContainer).getValue();

                    //note: we ignore all the Sonar suggestions here to keep the code compatible with Java 8
                    if (result instanceof DoubleResult) {
                        DoubleResult doubleResult = (DoubleResult) result;
                        vector[startIndex] = (float) doubleResult.doubleValue();
                        success = true;
                    } else if (result instanceof IntegerResult) {
                        IntegerResult integerResult = (IntegerResult) result;
                        vector[startIndex] = (float) integerResult.intValue();
                        success = true;
                    } else if (result instanceof DoubleArrayResult) {
                        DoubleArrayResult arrayResult = (DoubleArrayResult) result;
                        for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = (float) arrayResult.get(i);
                        }
                        success = true;
                    } else if (result instanceof IntegerArrayResult) {
                        IntegerArrayResult arrayResult = (IntegerArrayResult) result;
                        for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = (float) arrayResult.get(i);
                        }
                        success = true;
                    } else {
                        LOGGER.warn("Unexpected result type ", result.getClass().getName(), " for descriptor ", descriptor.name(), ". Expected DoubleResult, IntegerResult, DoubleArrayResult or IntegerArrayResult.");
                    }
                } else {
                    LOGGER.warn("Descriptor ", descriptor.name(), " not found in descriptorToCdkObjectMap. This should not happen.");
                }
            }
        } catch (InterruptedException exception) {
            LOGGER.warn("Interrupted while waiting for fingerprinter instance from pool. This should not happen.", exception);
            throw exception;
        } catch (Exception exception) {
            LOGGER.error("Descriptor ", descriptor.name(), " failed to calculate: ", exception.getMessage());
            return false;
        }
        return success;
    }

    /**
     * Initializes all fingerprint pools with the current {@link #fingerprintPoolSize}.
     * This method creates new BlockingQueues for each fingerprint type and populates them
     * with fingerprinter instances. Thread-safe for concurrent access.
     *
     * @return true if all pools were successfully initialized, false if any pool failed to initialize properly
     * (e.g., due to offer() failure, which should not happen under normal circumstances).
     */
    private static synchronized boolean initializeFingerprintPools() {
        boolean success = true;
        // Clear existing pools
        DescriptorCalculator.fingerprintPoolMap.clear();

        // PUBCHEM_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> pubchemPool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!pubchemPool.offer(new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()))){
                LOGGER.warn("Failed to add PubchemFingerprinter instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.PUBCHEM_FINGERPRINTER, pubchemPool);

        // CIRCULAR_FINGERPRINTER_ECFP_0 Pool
        BlockingQueue<IFingerprinter> ecfp0Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!ecfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP0, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter ECFP0 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0, ecfp0Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_0 Pool
        BlockingQueue<IFingerprinter> fcfp0Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!fcfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP0, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter FCFP0 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0, fcfp0Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_2 Pool
        BlockingQueue<IFingerprinter> ecfp2Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!ecfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP2, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter ECFP2 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2, ecfp2Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_2 Pool
        BlockingQueue<IFingerprinter> fcfp2Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!fcfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP2, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter FCFP2 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2, fcfp2Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_4 Pool
        BlockingQueue<IFingerprinter> ecfp4Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!ecfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter ECFP4 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, ecfp4Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_4 Pool
        BlockingQueue<IFingerprinter> fcfp4Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!fcfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP4, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter FCFP4 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4, fcfp4Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_6 Pool
        BlockingQueue<IFingerprinter> ecfp6Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!ecfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter ECFP6 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6, ecfp6Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_6 Pool
        BlockingQueue<IFingerprinter> fcfp6Pool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!fcfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6, DescriptorCalculator.circularFingerprintSize))) {
                LOGGER.warn("Failed to add CircularFingerprinter FCFP6 instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6, fcfp6Pool);

        // MACCS_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> maccsPool = new LinkedBlockingQueue<>(DescriptorCalculator.fingerprintPoolSize);
        for (int i = 0; i < DescriptorCalculator.fingerprintPoolSize; i++) {
            if (!maccsPool.offer(new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()))) {
                LOGGER.warn("Failed to add MACCSFingerprinter instance to pool. This should not happen.");
                success = false;
            }
        }
        DescriptorCalculator.fingerprintPoolMap.put(Descriptor.MACCS_FINGERPRINTER, maccsPool);

        return success;
    }

    /**
     * Sets fingerprint component values in vector.
     *
     * @param descriptor Fingerprint descriptor to calculate (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Vector of molecule (row in data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param startIndex Start index in vector to be filled with calculated components of descriptors
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * @throws InterruptedException if the current thread is interrupted while waiting to acquire a
     *                              fingerprinter instance from the pool; this is the root source of
     *                              {@link InterruptedException} in this class — {@link BlockingQueue#take()}
     *                              throws it if the thread is interrupted while blocked waiting for an
     *                              available fingerprinter, or if it is already interrupted when called
     */
    private static boolean calculateFingerprintFromPool(Descriptor descriptor, IAtomContainer atomContainer, float[] vector, int startIndex)
            throws InterruptedException {
        IFingerprinter fingerprinter = null;
        boolean success;
        try {
            BlockingQueue<IFingerprinter> pool = DescriptorCalculator.fingerprintPoolMap.get(descriptor);
            // take() blocks until a fingerprinter is available in the pool; throws InterruptedException
            // if the thread is interrupted while waiting (or already interrupted). InterruptedException is
            // neither CDKException nor RuntimeException, so it is not caught below and propagates naturally.
            fingerprinter = pool.take();

            IBitFingerprint bitFingerprint = fingerprinter.getBitFingerprint(atomContainer);
            for (int i = 0; i < descriptor.getDescriptorComponentNumber(); i++) {
                vector[startIndex + i] = bitFingerprint.get(i) ? 1.0f : 0.0f;
            }
            success = true;
        } catch (CDKException | RuntimeException exception) {
            LOGGER.warn("Failed to calculate: " + descriptor.name(), exception);
            return false;
        } finally {
            // Only return the fingerprinter if it was successfully taken from the pool.
            // If take() threw InterruptedException, fingerprinter is still null and offer(null)
            // would throw NullPointerException.
            if (fingerprinter != null && !DescriptorCalculator.fingerprintPoolMap.get(descriptor).offer(fingerprinter)){
                LOGGER.warn("Failed to return fingerprinter to pool: " + descriptor.name());
                success = false;
            }

        }
        return success;
    }

    /**
     * Sets the pool size for fingerprinter instances and reinitialized all fingerprint pools.
     * The default pool size is 4, which should be sufficient for regular users.
     * Increasing the pool size can improve performance in highly parallel environments
     * but will increase memory usage.
     * <p>
     * This method is thread-safe and will block until all pools are reinitialized.
     * Any fingerprinter instances currently in use will be returned to the old pools
     * and will eventually be garbage collected.
     *
     * @param aPoolSize The new pool size for fingerprinter instances (must be greater than 0)
     * @throws IllegalArgumentException if aPoolSize is less than or equal to 0
     */
    public static synchronized boolean setFingerprintPoolSize(int aPoolSize) throws IllegalArgumentException {
        if (aPoolSize <= 0) {
            throw new IllegalArgumentException("Descriptor.setFingerprintPoolSize: aPoolSize must be greater than 0.");
        }
        try {
            DescriptorCalculator.fingerprintPoolSize = aPoolSize;
            DescriptorCalculator.initializeFingerprintPools();
            return true;
        } catch (Exception exception){
            LOGGER.error("Failed to set fingerprint pool size and reinitialize pools: ", exception.getMessage());
            return false;
        }

    }

    /**
     * Returns the current pool size for fingerprinter instances.
     *
     * @return Current fingerprint pool size
     */
    public static synchronized int getFingerprintPoolSize() {
        return DescriptorCalculator.fingerprintPoolSize;
    }

    /**
     * Sets the size for circular fingerprints and reinitializes the corresponding pools.
     * This method is thread-safe and will block until all pools are reinitialized.
     * Any fingerprinter instances currently in use will be returned to the old pools
     * and will eventually be garbage collected.
     *
     * @param aSize The new size for circular fingerprints (must be greater than 0)
     * @return true if successful, false otherwise
     * @throws IllegalArgumentException if aSize is less than or equal to 0
     */
    public static synchronized boolean setCircularFingerprintSize(int aSize) throws IllegalArgumentException {
        if (aSize <= 0) {
            throw new IllegalArgumentException("Descriptor.setCircularFingerprintSize: aSize must be greater than 0.");
        }
        try {
            DescriptorCalculator.circularFingerprintSize = aSize;
            Descriptor.setCircularFingerprintComponentNumber(aSize);
            DescriptorCalculator.initializeFingerprintPools();
            return true;
        } catch (Exception exception) {
            LOGGER.error("Failed to set circular fingerprint size and reinitialize pools: ", exception.getMessage());
            return false;
        }
    }

    /**
     * Returns the current size used for circular fingerprints.
     *
     * @return Current circular fingerprint size
     */
    public static synchronized int getCircularFingerprintSize() {
        return DescriptorCalculator.circularFingerprintSize;
    }


    /**
     * Uses a specified aromaticity model to modify molecule.
     * The method performs a complete workflow of:
     * <ol>
     *     <li>Suppresses explicit hydrogens.</li>
     *     <li>Perceives atom types and configures atoms.</li>
     *     <li>Clears existing aromaticity flags.</li>
     *     <li>Detects rings.</li>
     *     <li>Applies the specified {@link ElectronDonation} model.</li>
     * </ol>
     * Note: This method changes the input molecule by applying the specified aromaticity model and reconfiguring its hydrogens.
     *
     * @param molecule Molecule that will be modified (IS CHANGED)
     * @param electronDonationModel The aromaticity model that will be used for aromaticity detection
     * @throws CDKException if adding implicit hydrogen atoms or detection of atom types fails
     */
    public static void setAromaticity(
            IAtomContainer molecule,
            ElectronDonation electronDonationModel
    ) throws CDKException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (molecule.isEmpty()) {
            return;
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("Aromaticity model must not be null");
        }
        AtomContainerManipulator.normalizeHydrogens(molecule, HydrogenState.Minimal);
        // Needed for VABC Descriptor and CDK_AtomTypes
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
        // Clears all aromatic flags before applying the aromaticity model.
        Aromaticity.clear(molecule);
        Cycles.markRingAtomsAndBonds(molecule);
        Aromaticity.apply(electronDonationModel, molecule);
    }

    /**
     * Calculates a single descriptor for the given molecule and returns the result as a float array.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for
     * single descriptor calculations. If the calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     * <p>
     * <b>Important:</b> The molecule must have aromaticity already perceived if required by the descriptor.
     * Use {@link #setAromaticity(IAtomContainer, ElectronDonation)} before calling this method if needed.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param molecule The molecule to calculate the descriptor for (must not be null);
     *                  aromaticity must be perceived beforehand if required; if it is empty, an empty array is returned
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link Descriptor#getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor, IAtomContainer molecule)  {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }
        if (descriptor.requires3DCoordinates() && !GeometryUtil.has3DCoordinates(molecule)) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("DescriptorCalculator.calculateDescriptor: Descriptor " + descriptor.getName() + " requires 3D coordinates, but the molecule does not have 3D coordinates.");
            return result;
        }

        float[] result = new float[descriptor.getDescriptorComponentNumber()];
        try {
            boolean success;
            if (descriptor.needsExplicitHydrogens()) {
                IAtomContainer moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(molecule);
                success = DescriptorCalculator.calculate(descriptor, moleculeWithExplicitHydrogens, result, 0);
            } else {
                success = DescriptorCalculator.calculate(descriptor, molecule, result, 0);
            }
            if (!success) {
                Arrays.fill(result, Float.NaN);
            }
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("Interrupted while calculating descriptor " + descriptor.getName(), exception);
        } catch (Exception exception) {
            // Fill with NaN on failure and log the error
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
        }
        return result;
    }

    /**
     * Calculates a single descriptor for the given molecule with automatic aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically perceives aromaticity using the specified electron donation model before calculating the descriptor.
     * If aromaticity perception or calculation fails, all values in the result array will be set to {@link Float#NaN}.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param molecule The molecule to calculate the descriptor for (must not be null); if it is empty, an empty array is returned
     * @param electronDonationModel The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link Descriptor#getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor,
                                              IAtomContainer molecule,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: electronDonationModel must not be null.");
        }

        try {
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
            return DescriptorCalculator.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
            return result;
        }
    }

    /**
     * Calculates a single descriptor for the given SMILES string with automatic parsing and aromaticity perception.
     * <p>
     * This method is a convenience wrapper around the internal {@link #calculate} method for single descriptor calculations.
     * The method automatically parses the SMILES string into a molecule,
     * perceives aromaticity using the specified electron donation model, and then calculates the descriptor.
     * If parsing, aromaticity perception, or calculation fails, all values in the result array
     * will be set to {@link Float#NaN}.
     *
     * @param descriptor The descriptor to calculate (must not be null)
     * @param smilesString The SMILES string representing the molecule (must not be null); if it is blank, an empty array is returned
     * @param electronDonationModel The electron donation model to use for aromaticity perception (must not be null)
     * @return A float array containing the calculated descriptor components. Length equals
     *         {@link Descriptor#getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor,
                                              String smilesString,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: descriptor must not be null.");
        }
        if (smilesString == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: smilesString must not be null.");
        }
        if (smilesString.trim().isEmpty()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("DescriptorCalculator.calculateDescriptor: electronDonationModel must not be null.");
        }
        if (descriptor.requires3DCoordinates()) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn("DescriptorCalculator.calculateDescriptor: Descriptor " + descriptor.getName() + " requires 3D coordinates, which cannot be computed from a SMILES string.");
            return result;
        }

        try {
            IAtomContainer molecule = DescriptorCalculator.SMILES_PARSER.parseSmiles(smilesString);
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
            return DescriptorCalculator.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            DescriptorCalculator.LOGGER.warn(String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
            return result;
        }
    }


    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecule batches. If parallel computation is used, the atom
     * container batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via {@link #setFingerprintPoolSize(int)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainerArray Array of molecules. Note: atomContainerArray[i] corresponds to matrix[i] data
     *                             vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                atomContainerArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param batchSize Number of molecules to process in each batch
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (isParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN or a global exception occurred
     * @throws IllegalArgumentException if the matrix dimensions are invalid or if batch size is &lt;= 0
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected batch is aborted silently
     */
    public static boolean setDescriptorsForMoleculesByBatchParallelization(
            Descriptor[] descriptors,
            IAtomContainer[] atomContainerArray,
            float[][] matrix,
            int startIndex,
            int batchSize,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesByBatchParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return true;
        }
        if (!DescriptorCalculator.validateAtomContainerArray(atomContainerArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given atom container array is empty, calculation aborted.");
            return true;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                if (atomContainerArray[0] == null || !GeometryUtil.has3DCoordinates(atomContainerArray[0])) {
                    DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, but the first molecule does not have 3D coordinates. Calculation aborted.");
                    return true;
                }
                break;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        DescriptorCalculator.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

        int numberOfMolecules = atomContainerArray.length;
        int numberOfBatches = (int) Math.ceil((double) numberOfMolecules / batchSize);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += descriptors[i].getDescriptorComponentNumber();
        }

        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                IntStream.range(0, numberOfBatches).parallel().forEach(batchIndex -> {
                    int batchStart = batchIndex * batchSize;
                    int batchEnd = Math.min(batchStart + batchSize, numberOfMolecules);

                    for (int i = batchStart; i < batchEnd; i++) {
                        try {
                            boolean success = DescriptorCalculator.setDescriptorsForSingleMolecule(
                                    descriptors,
                                    atomContainerArray[i],
                                    matrix[i],
                                    startIndices,
                                    i,
                                    nanPositionsList
                            );
                            if (!success) {
                                hasNaN.set(true);
                            }
                        } catch (InterruptedException exception) {
                            Thread.currentThread().interrupt(); // Preserve interrupt status
                            return;
                        } catch (Exception exception) {
                            hasNaN.set(true);
                            DescriptorCalculator.LOGGER.warn(
                                    String.format("Descriptor.setDescriptorsForMoleculesByBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!DescriptorCalculator.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, nanPositionsList)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("Descriptor.setDescriptorsForMoleculesByBatchParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "Descriptor.setDescriptorsForMoleculesByBatchParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecule batches. If parallel computation is used, the smiles string
     * batches (molecules) are distributed onto parallel thread, one for each molecule batch, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via {@link #setFingerprintPoolSize(int)}.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesStringArray Array of molecule smiles strings. Note: moleculeSmilesStringArray[i] corresponds to matrix[i] data
     *                                   vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                moleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param batchSize Number of molecules to process in each batch
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                                Aromaticity.Model.Daylight is used as default.
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (isParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN or a global exception occurred
     * @throws IllegalArgumentException if the matrix dimensions are invalid or if batch size is &lt;= 0
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected batch is aborted silently
     */
    public static boolean setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
            Descriptor[] descriptors,
            String[] moleculeSmilesStringArray,
            float[][] matrix,
            int startIndex,
            int batchSize,
            ElectronDonation electronDonationModel,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculeBySmilesStringsBatchParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return true;
        }
        if (!DescriptorCalculator.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given SMILES string array is empty, calculation aborted.");
            return true;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, which cannot be computed from SMILES strings. Calculation aborted.");
                return true;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        DescriptorCalculator.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

        int numberOfMolecules = moleculeSmilesStringArray.length;
        int numberOfBatches = (int) Math.ceil((double) numberOfMolecules / batchSize);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += descriptors[i].getDescriptorComponentNumber();
        }

        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                IntStream.range(0, numberOfBatches).parallel().forEach(batchIndex -> {
                    int batchStart = batchIndex * batchSize;
                    int batchEnd = Math.min(batchStart + batchSize, numberOfMolecules);

                    for (int i = batchStart; i < batchEnd; i++) {
                        try {
                            boolean success = DescriptorCalculator.setDescriptorsForSingleMoleculeSmilesString(
                                    descriptors,
                                    moleculeSmilesStringArray[i],
                                    matrix[i],
                                    startIndices,
                                    i,
                                    electronDonationModel,
                                    nanPositionsList
                            );
                            if (!success) {
                                hasNaN.set(true);
                            }
                        } catch (InterruptedException exception) {
                            Thread.currentThread().interrupt(); // Preserve interrupt status
                            return;
                        } catch (Exception exception) {
                            hasNaN.set(true);
                            DescriptorCalculator.LOGGER.warn(
                                    String.format("Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!DescriptorCalculator.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, nanPositionsList)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via {@link #setFingerprintPoolSize(int)}
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainerArray Array of molecules. Note: atomContainerArray[i] corresponds to matrix[i] data
     *                             vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                atomContainerArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (isParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN or a global exception occurred
     * @throws IllegalArgumentException if the matrix dimensions are invalid
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected molecule is skipped silently
     */
    public static boolean setDescriptorsForMoleculesByMoleculeParallelization(
            Descriptor[] descriptors,
            IAtomContainer[] atomContainerArray,
            float[][] matrix,
            int startIndex,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesByMoleculeParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return true;
        }
        if (!DescriptorCalculator.validateAtomContainerArray(atomContainerArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given atom container array is empty, calculation aborted.");
            return true;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                if (atomContainerArray[0] == null || !GeometryUtil.has3DCoordinates(atomContainerArray[0])) {
                    DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, but the first molecule does not have 3D coordinates. Calculation aborted.");
                    return true;
                }
                break;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        DescriptorCalculator.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += descriptors[i].getDescriptorComponentNumber();
        }

        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                // Advise by Oracle: Parallel streams should use the common Fork-join pool
                IntStream.range(0, atomContainerArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                boolean success = DescriptorCalculator.setDescriptorsForSingleMolecule(
                                        descriptors,
                                        atomContainerArray[i],
                                        matrix[i],
                                        startIndices,
                                        i,
                                        nanPositionsList
                                );
                                if (!success) {
                                    hasNaN.set(true);
                                }
                            } catch (InterruptedException exception) {
                                Thread.currentThread().interrupt(); // Preserve interrupt status
                            } catch (Exception exception) {
                                hasNaN.set(true);
                                DescriptorCalculator.LOGGER.warn(
                                        String.format("Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: One descriptor calculation caused an exception, molecule index: %d.", i),
                                        exception
                                );
                            }
                        }
                );
            } else {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    try {
                        if (!DescriptorCalculator.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, nanPositionsList)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: Global exception occurred in descriptor calculation: ",
                    exception
            );
            return false;
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecules. If parallel computation is used, the smiles strings
     * (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via {@link #setFingerprintPoolSize(int)}
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesStringArray Array of molecule smiles strings. Note: moleculeSmilesStringArray[i] corresponds to matrix[i] data
     *                                   vector, i.e. the molecules define the rows of the matrix (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                moleculeSmilesStringArray[i]. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     *                                Aromaticity.Model.Daylight is used as default.
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (isParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN or a global exception occurred
     * @throws IllegalArgumentException if the matrix dimensions are invalid
     * @throws InterruptedException if the current thread is interrupted during sequential fingerprint
     *                              calculation; in parallel mode, interruption is handled internally
     *                              by restoring the thread's interrupt flag via
     *                              {@link Thread#interrupt()} and the affected molecule is skipped silently
     */
    public static boolean setDescriptorsForMoleculesBySmilesStringParallelization(
            Descriptor[] descriptors,
            String[] moleculeSmilesStringArray,
            float[][] matrix,
            int startIndex,
            ElectronDonation electronDonationModel,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException, InterruptedException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesBySmilesStringParallelization";
        if (!DescriptorCalculator.validateDescriptors(descriptors, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return true;
        }
        if (!DescriptorCalculator.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            DescriptorCalculator.LOGGER.warn(methodName + " : Given SMILES string array is empty, calculation aborted.");
            return true;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                DescriptorCalculator.LOGGER.warn(methodName + " : At least one descriptor requires 3D coordinates, which cannot be computed from SMILES strings. Calculation aborted.");
                return true;
            }
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        DescriptorCalculator.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

        int numberOfMolecules = moleculeSmilesStringArray.length;
        int[] startIndices = new int[descriptors.length];
        for (int i = 0; i < descriptors.length; i++) {
            startIndices[i] = startIndex;
            startIndex += descriptors[i].getDescriptorComponentNumber();
        }
        AtomicBoolean hasNaN = new AtomicBoolean(false);
        try {
            if (isParallelCalculation) {

                // Advise by Oracle: Parallel streams should use the common Fork-join pool
                IntStream.range(0, moleculeSmilesStringArray.length).parallel().forEach(
                        i ->
                        {
                            try {
                                boolean success = DescriptorCalculator.setDescriptorsForSingleMoleculeSmilesString(
                                        descriptors,
                                        moleculeSmilesStringArray[i],
                                        matrix[i],
                                        startIndices,
                                        i,
                                        electronDonationModel,
                                        nanPositionsList
                                );
                                if (!success) {
                                    hasNaN.set(true);
                                }
                            } catch (InterruptedException exception) {
                                Thread.currentThread().interrupt(); // Preserve interrupt status
                            } catch (Exception exception) {
                                hasNaN.set(true);
                                DescriptorCalculator.LOGGER.warn(
                                        String.format("Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Exception in molecule index: %d", i),
                                        exception
                                );
                            }

                        }
                );
            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    try {
                        if (!DescriptorCalculator.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, nanPositionsList)) {
                            hasNaN.set(true);
                        }
                    } catch (InterruptedException exception) {
                        throw exception;
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        DescriptorCalculator.LOGGER.warn(
                                String.format("Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            DescriptorCalculator.LOGGER.warn(
                    "Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Global exception occurred in descriptor calculation.",
                    exception
            );
            return false;
        }
        return !hasNaN.get();
    }

    /**
     * Sets calculated descriptor components in vector (that corresponds to atomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in vector to be filled with calculated components of descriptors (each
     *                      descriptor in descriptors has its dedicated start index here)
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws CloneNotSupportedException if copying the molecule fails
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #setDescriptor})
     */
    private static boolean setDescriptorsForSingleMolecule (
            Descriptor[] descriptors,
            IAtomContainer atomContainer,
            float[] vector,
            int[] aStartIndices,
            int moleculeIndex,
            List<int[]> nanPositionsList
    ) throws CloneNotSupportedException, InterruptedException {
        IAtomContainer moleculeWithExplicitHydrogens = null;
        for (Descriptor descriptor : descriptors) {
            if (descriptor.needsExplicitHydrogens()){
                moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(atomContainer);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : atomContainer;
            if (!DescriptorCalculator.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
                isSuccessful = false;
            }
        }
        return isSuccessful;
    }

    /**
     * Sets calculated descriptor components in vector (that corresponds to atomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param moleculeSmilesString Molecule SMILES String (IS NOT CHANGED)
     * @param vector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in vector to be filled with calculated components of descriptors (each
     *                      descriptor in descriptors has its dedicated start index here)
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @param electronDonationModel An Aromaticity model that is applied to every molecule. NOTE: Can be null, then
     * Aromaticity.Model.Daylight is used as default.
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws CloneNotSupportedException if copying the molecule fails
     * @throws CDKException if SMILES parsing or aromaticity detection fails
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #setDescriptor})
     */
    private static boolean setDescriptorsForSingleMoleculeSmilesString (
            Descriptor[] descriptors,
            String moleculeSmilesString,
            float[] vector,
            int[] aStartIndices,
            int moleculeIndex,
            ElectronDonation electronDonationModel,
            List<int[]> nanPositionsList
    ) throws CDKException, CloneNotSupportedException, InterruptedException {
        IAtomContainer molecule = DescriptorCalculator.SMILES_PARSER.parseSmiles(moleculeSmilesString);
        if (electronDonationModel == null) {
            DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
        } else {
            DescriptorCalculator.setAromaticity(molecule, electronDonationModel);
        }

        IAtomContainer moleculeWithExplicitHydrogens = null;
        for (Descriptor descriptor : descriptors) {
            if (descriptor.needsExplicitHydrogens()){
                moleculeWithExplicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(molecule);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : molecule;
            if (!DescriptorCalculator.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
                isSuccessful = false;
            }
        }
        return isSuccessful;
    }

    /**
     * Sets component values of descriptor for atomContainer in vector beginning with startIndex.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Vector of molecule (row of data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param startIndex Start index in vector to be filled with calculated components of descriptors
     * @param moleculeIndex Index of the current molecule being processed (row index in data matrix)
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws InterruptedException if the current thread is interrupted during fingerprint calculation
     *                              (propagated from {@link #calculate}); all other exceptions are caught,
     *                              logged, and converted to NaN values in vector
     */
    private static boolean setDescriptor(
            Descriptor descriptor,
            IAtomContainer atomContainer,
            float[] vector,
            int startIndex,
            int moleculeIndex,
            List<int[]> nanPositionsList
    ) throws InterruptedException {
        try {
            boolean success = DescriptorCalculator.calculate(descriptor, atomContainer, vector, startIndex);
            if (!success) {
                int numComponents = descriptor.getDescriptorComponentNumber();
                for (int i = 0; i < numComponents; i++) {
                    vector[startIndex + i] = Float.NaN;
                    if (nanPositionsList != null) {
                        nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                    }
                }
                return false;
            }

            // Check for NaN values in the calculated result and track them
            int numComponents = descriptor.getDescriptorComponentNumber();
            return !DescriptorCalculator.checkAndTrackNaNValues(vector, startIndex, numComponents, moleculeIndex, nanPositionsList);
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            int numComponents = descriptor.getDescriptorComponentNumber();
            for (int i = 0; i < numComponents; i++) {
                vector[startIndex + i] = Float.NaN;
                // Track NaN position if nanPositionsList is provided
                if (nanPositionsList != null) {
                    nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                }
            }
            DescriptorCalculator.LOGGER.warn(
                    String.format("Descriptor.setDescriptor: An exception occurred while calculating descriptor %s for molecule index %d.",
                            descriptor,
                            moleculeIndex),
                    exception
            );
            return false;
        }
    }

    /**
     * Helper method to check for NaN values in a calculated descriptor result and track their positions.
     * This method is thread-safe when used with a synchronized LinkedList. Note that not the entire data vector is
     * checked but only the positions [startIndex, startIndex + numComponents].
     *
     * @param vector The vector containing calculated descriptor values
     * @param startIndex The start index in the vector for this descriptor
     * @param numComponents The number of components for this descriptor
     * @param moleculeIndex The index of the current molecule
     * @param nanPositionsList List to track NaN positions (can be null if NaN positions should not be tracked);
     *                      must be thread-safe for parallel access; will be filled with int[] {[moleculeIndex, componentIndex]}
     *                      pairs of NaN value positions (i.e. row and column index in the data matrix)
     * @return true if any NaN values were found, false otherwise
     */
    static boolean checkAndTrackNaNValues(
            float[] vector,
            int startIndex,
            int numComponents,
            int moleculeIndex,
            List<int[]> nanPositionsList
    ) {
        boolean foundNaN = false;
        for (int i = 0; i < numComponents; i++) {
            if (Float.isNaN(vector[startIndex + i])) {
                foundNaN = true;
                // Track NaN position if nanPositionsList is provided
                // Note: nanPositionsList must be thread-safe for parallel access
                if (nanPositionsList != null) {
                    // Synchronize the add operation to ensure thread safety
                    nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                }
            }
        }
        return foundNaN;
    }

    /**
     * Validates descriptor array input. Throws NullPointerExceptions if the array or one of its elements is null.
     * Returns false if the array is empty.
     *
     * @param descriptorArray the descriptor array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if descriptor array or one of its elements is null
     * @return false if descriptor array is empty; true otherwise
     */
    static boolean validateDescriptors(Descriptor[] descriptorArray, String methodName) {
        if (descriptorArray == null) {
            throw new NullPointerException(methodName + ": Given descriptor array is null");
        }
        if (descriptorArray.length == 0) {
            return false;
        }
        for (Descriptor descriptor : descriptorArray) {
            if (descriptor == null) {
                throw new NullPointerException(methodName + ": A single descriptor in descriptors is null.");
            }
        }
        return true;
    }

    /**
     * Validates an atom container array input and its contents. Throws NullPointerExceptions if the array or one of
     * its elements is null. Returns false if the array is empty. Note that the molecules in the array are allowed to be empty.
     *
     * @param atomContainerArray the atom container array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if molecule array or one of its elements is null
     * @return false if the atom container array is empty; true otherwise
     */
    static boolean validateAtomContainerArray(IAtomContainer[] atomContainerArray, String methodName) {
        if (atomContainerArray == null) {
            throw new NullPointerException(methodName + ": atomContainerArray is null.");
        }
        if (atomContainerArray.length == 0) {
            return false;
        }
        for (IAtomContainer molecule : atomContainerArray) {
            if (molecule == null) {
                throw new NullPointerException(methodName + ": A single molecule in atomContainerArray is null.");
            }
            //do nothing if a molecule is empty
        }
        return true;
    }

    /**
     * Validates a SMILES string array input and its content. Throws NullPointerExceptions if the array or one of
     * its elements is null. Returns false if the array is empty.
     *
     * @param moleculeSmilesStringArray the molecule SMILES string array to validate
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if the SMILES string array or one of its elements is null
     * @return false if the SMILES string array is empty; true otherwise
     */
    private static boolean validateSmilesStringArray(String[] moleculeSmilesStringArray, String methodName) {
        if (moleculeSmilesStringArray == null) {
            throw new NullPointerException(methodName + ": moleculeSmilesStringArray is null.");
        }
        if (moleculeSmilesStringArray.length == 0) {
            return false;
        }
        for (String molecule : moleculeSmilesStringArray) {
            if (molecule == null) {
                throw new NullPointerException(methodName + ": A single SMILES string in moleculeSmilesStringArray is null.");
            }
            //do nothing if a String is empty
        }
        return true;
    }

    /**
     * Validates matrix dimensions and start index. Throws a NullPointerException if the matrix or one of its rows is
     * null. Returns false if the matrix is empty. Note that some parameters like {@code descriptors} have their own validation methods.
     *
     * @param matrix the matrix to validate
     * @param descriptors the descriptors intended to be calculated for validating the matrix dimensions (min. nr. of columns)
     * @param moleculeArray the input molecule array (either an {@code IAtomContainer[]} or a {@code String[]}) for
     *                       validating the matrix dimensions (exact(!) nr. of rows)
     * @param startIndex the start index in the matrix (starting there, enough columns must be left to fit all descriptor results)
     * @param methodName the calling method name for error messages
     * @throws NullPointerException if the matrix or one of its rows is null
     * @throws IllegalArgumentException if matrix dimensions are invalid
     */
    static void validateMatrix(
            float[][] matrix,
            Descriptor[] descriptors,
            Object[] moleculeArray,
            int startIndex,
            String methodName
    ) {
        if (matrix == null) {
            throw new NullPointerException(methodName + ": matrix is null.");
        }
        if (matrix.length == 0) {
            throw new IllegalArgumentException(methodName + ": matrix is of length 0.");
        }
        if (matrix.length != moleculeArray.length) {
            throw new IllegalArgumentException(methodName + ": matrix and moleculeArray must have the same length.");
        }
        if (!(moleculeArray instanceof IAtomContainer[]) && !(moleculeArray instanceof String[])) {
            throw new IllegalArgumentException(methodName + ": moleculeArray must be an IAtomContainer[] or a String[].");
        }
        int colNum = -1;
        for (float[] vector : matrix) {
            if (vector == null) {
                throw new NullPointerException(methodName + ": A vector in matrix is null.");
            }
            if (vector.length == 0) {
                throw new IllegalArgumentException(methodName + ": A vector in matrix has length 0.");
            }
            if (colNum == -1) {
                colNum = vector.length;
            }
            if (vector.length != colNum) {
                throw new IllegalArgumentException(methodName + ": All vectors in matrix must have the same length.");
            }
            if (startIndex > vector.length) {
                throw new IllegalArgumentException(methodName + ": startIndex is greater than or equal to vector length.");
            }
            int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
            if (startIndex + numberOfComponents > vector.length) {
                throw new IllegalArgumentException(methodName + ": Not enough space in vector for descriptors.");
            }
        }
    }

    /**
     * Creates a deep copy of the input molecule.
     * <p>
     * Notes:
     * <ul>
     *     <li>The copy is built with the same {@link org.openscience.cdk.interfaces.IChemObjectBuilder}
     *         as the input, using the in-container factory methods
     *         {@code newAtom(int, int)} / {@code newBond(IAtom, IAtom, IBond.Order)} which are
     *         faster than the dynamic {@code newInstance(...)} dispatch.</li>
     *     <li>Implicit hydrogen counts are preserved exactly (including a possible {@code null}).</li>
     *     <li>Bond endpoints are resolved through an explicit original-to-copy atom map, so no
     *         assumption is made about identical atom indices in both containers.</li>
     *     <li>If needed, atom types must be perceived and configured manually after creation.</li>
     * </ul>
     *
     * @param molecule Source molecule to be copied (NOT MODIFIED)
     * @return New independent instance of the molecule
     * @throws NullPointerException     if {@code molecule} is {@code null}
     * @throws CloneNotSupportedException if the molecule cannot be properly copied
     */
    static IAtomContainer copyMolecule(IAtomContainer molecule)
            throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        try {
            IAtomContainer moleculeCopy = molecule.getBuilder().newInstance(IAtomContainer.class);
            // Molecule-level properties and identifier
            if (molecule.getProperties() != null) {
                moleculeCopy.addProperties(new HashMap<>(molecule.getProperties()));
            }
            // Explicit original -> copy maps; used for bonds, stereo elements, electrons and lone pairs.
            int atomCapacity = Math.max(16, (int) (molecule.getAtomCount() / 0.75f) + 1);
            int bondCapacity = Math.max(16, (int) (molecule.getBondCount() / 0.75f) + 1);
            Map<IAtom, IAtom> atomMap = new HashMap<>(atomCapacity);
            Map<IBond, IBond> bondMap = new HashMap<>(bondCapacity);

            // ---- Atoms ----
            for (IAtom atom : molecule.atoms()) {
                Integer atomicNumber = atom.getAtomicNumber();
                Integer implicitHCount = atom.getImplicitHydrogenCount();
                // fast in-container creation (avoids reflective newInstance dispatch)
                IAtom newAtom = moleculeCopy.newAtom(
                        atomicNumber != null ? atomicNumber : 0,          // 0 == wildcard element
                        implicitHCount != null ? implicitHCount : 0);
                newAtom.setSymbol(atom.getSymbol());
                newAtom.setImplicitHydrogenCount(implicitHCount);         // preserves a possible null
                newAtom.setMassNumber(atom.getMassNumber());
                newAtom.setFormalCharge(atom.getFormalCharge());
                newAtom.setCharge(atom.getCharge());
                newAtom.setValency(atom.getValency());
                newAtom.setHybridization(atom.getHybridization());
                newAtom.setAtomTypeName(atom.getAtomTypeName());
                // Coordinates (needed for stereo perception / preservation)
                Point2d p2d = atom.getPoint2d();
                if (p2d != null) {
                    newAtom.setPoint2d(new Point2d(p2d));
                }
                Point3d p3d = atom.getPoint3d();
                if (p3d != null) {
                    newAtom.setPoint3d(new Point3d(p3d));
                }
                Point3d fp3d = atom.getFractionalPoint3d();
                if (fp3d != null) {
                    newAtom.setFractionalPoint3d(new Point3d(fp3d));
                }
                // Flags
                newAtom.setIsAromatic(atom.isAromatic());
                newAtom.setIsInRing(atom.isInRing());
                // Full generic property map (carries the SRU unique atom index + detection markers)
                if (atom.getProperties() != null) {
                    newAtom.addProperties(new HashMap<>(atom.getProperties()));
                }
                atomMap.put(atom, newAtom);
            }

            // ---- Bonds ----
            for (IBond bond : molecule.bonds()) {
                IAtom beginCopy = atomMap.get(bond.getBegin());
                IAtom endCopy = atomMap.get(bond.getEnd());
                IBond newBond = moleculeCopy.newBond(beginCopy, endCopy, bond.getOrder());
                newBond.setDisplay(bond.getDisplay());
                newBond.setIsAromatic(bond.isAromatic());
                newBond.setIsInRing(bond.isInRing());
                if (bond.getProperties() != null) {
                    newBond.addProperties(new HashMap<>(bond.getProperties()));
                }
                bondMap.put(bond, newBond);
            }

            // ---- Single electrons ----
            for (ISingleElectron singleElectron : molecule.singleElectrons()) {
                ISingleElectron newSingleElectron = molecule.getBuilder().newInstance(
                        ISingleElectron.class, atomMap.get(singleElectron.getAtom()));
                if (singleElectron.getProperties() != null) {
                    newSingleElectron.addProperties(new HashMap<>(singleElectron.getProperties()));
                }
                moleculeCopy.addSingleElectron(newSingleElectron);
            }

            // ---- Stereo elements (preserve stereochemistry, re-mapped onto the copy) ----
            for (IStereoElement<?, ?> stereoElement : molecule.stereoElements()) {
                moleculeCopy.addStereoElement(stereoElement.map(atomMap, bondMap));
            }

            return moleculeCopy;
        } catch (Exception exception) {
            throw new CloneNotSupportedException("Could not clone molecule: " + exception.getMessage());
        }
    }

    /**
     * Creates a molecule copy with all implicit hydrogen atoms made explicit.
     *
     * @param molecule Source molecule with implicit hydrogens (NOT MODIFIED)
     * @return New molecule (copy of the parameter) with all hydrogens made explicit
     * @throws CloneNotSupportedException If the molecule cannot be properly processed
     */
    static IAtomContainer createMoleculeWithExplicitHydrogens(
            IAtomContainer molecule
    ) throws NullPointerException, IllegalArgumentException, CloneNotSupportedException {
        // Checks
        if (molecule == null) {
            throw new NullPointerException("Input molecule must not be null");
        }
        if (molecule.isEmpty()) {
            //returns empty atom container
            return DescriptorCalculator.copyMolecule(molecule);
        }
        try {
            // Create a deep copy of the molecule first
            IAtomContainer moleculeCopy = DescriptorCalculator.copyMolecule(molecule);
            // Add explicit hydrogen atoms
            AtomContainerManipulator.normalizeHydrogens(moleculeCopy, HydrogenState.Explicit);
            return moleculeCopy;
        } catch (Exception exception) {
            throw new CloneNotSupportedException("Could not create molecule with explicit hydrogens: " + exception.getMessage());
        }
    }
}
