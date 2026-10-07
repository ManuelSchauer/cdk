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

package de.unijena.cheminf.clustering.desccalc;

import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.ISingleElectron;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.qsar.IMolecularDescriptor;
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
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.HydrogenState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.IntStream;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;

//TODO: add a note on which descriptors are included here (all IMolecularDescriptor implementing classes except for
// those that require 3D coordinates, correct? Actually, we could still add them along with a new field "requires3DCoordinates")
/**
 * Descriptor related calculations based on the CDK for the enrichment of data vectors.
 * <p>
 *     Each enum constant carries six fields that describe its characteristics and can be queried via the
 *     corresponding getter methods:
 *     <ul>
 *         <li><b>isFast</b> ({@code boolean}) – whether this descriptor can be calculated quickly.
 *             The classification is based on a performance benchmark: descriptors that complete computation
 *             for 10,000 molecules in under one second are considered fast.
 *             Slow descriptors (e.g., {@link #CHI_CHAIN}, {@link #CHI_CLUSTER}, {@link #CHI_PATH_CLUSTER},
 *             {@link #CHI_PATH}, {@link #WEIGHTED_PATH}) require more intensive computation.
 *             Use {@link #isFast()} to query this field.</li>
 *         <li><b>isSafe</b> ({@code boolean}) – whether this descriptor reliably produces a valid result
 *             for a valid molecule input. Unsafe descriptors (e.g., {@link #BCUT}, {@link #KAPPA_SHAPE_INDICES})
 *             may produce {@code NaN} values for certain molecular inputs they are not developed for.
 *             Use {@link #isSafe()} to query this field.</li>
 *         <li><b>isFingerprint</b> ({@code boolean}) – whether this descriptor generates a binary
 *             fingerprint bit vector rather than a set of numerical descriptor values. Fingerprint
 *             descriptors (e.g., {@link #PUBCHEM_FINGERPRINTER}, {@link #MACCS_FINGERPRINTER},
 *             all {@code CIRCULAR_FINGERPRINTER_*} variants) are handled via a thread-safe pool of
 *             {@link IFingerprinter} instances internally.
 *             Use {@link #isFingerprint()} to query this field.</li>
 *         <li><b>needsExplicitHydrogens</b> ({@code boolean}) – whether this descriptor requires
 *             explicit hydrogen atoms to be present in the molecule before calculation.
 *             Examples include {@link #A_LOG_P}, {@link #X_LOG_P}, {@link #KAPPA_SHAPE_INDICES},
 *             and {@link #HYBRIDIZATION_RATIO}.
 *             Use {@link #needsExplicitHydrogens()} to query this field.</li>
 *         <li><b>descriptorComponentNumber</b> ({@code int}) – the number of individual numerical
 *             values produced by this descriptor. Most descriptors return a single value (e.g.,
 *             {@link #MOLECULAR_WEIGHT}, {@link #TPSA}), while multi-value descriptors return arrays
 *             (e.g., {@link #BCUT} returns 6, {@link #KAPPA_SHAPE_INDICES} returns 3,
 *             {@link #CHI_PATH} returns 16, {@link #KIER_HALL_SMARTS} returns 79,
 *             {@link #PUBCHEM_FINGERPRINTER} returns 881).
 *             Use {@link #getDescriptorComponentNumber()} to query this field.</li>
 *         <li><b>name</b> ({@code String}) – a human-readable display name for this descriptor,
 *             intended for labeling output columns or user-facing representations (e.g.,
 *             {@code "Molecular Weight"}, {@code "TPSA"}, {@code "ALogP"}).
 *             Use {@link #getName()} to query this field.</li>
 *     </ul>
 * <p>
 *     TODO: correct this to 4 after removing the "new" method<br>
 *     There are 5 different "public static boolean setDescriptorsForMolecules...()" methods with different forms of
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
 *     Example usage of the Descriptor class:
 * </p>
 * <pre>{@code
 * // Preprocessing
 * // A: Molecule preprocessing:
 * // 1: Create a molecule
 * IAtomContainer molecule = ...;
 *
 * // 2. IMPORTANT: Perform aromaticity perception before descriptor calculation
 * Descriptor.setAromaticity(molecule);
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
 * Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(); // e.g., user-defined selection
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
 * boolean result = Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
 *                      descriptors,
 *                      molecules,
 *                      matrix,
 *                      startIndex,
 *                      isParallelCalculation,
 *                      nanPositionsListSynchronized
 * );
 * // you can also use Descriptor.setDescriptorsForMoleculesByBatchParallelization(...) -> batchSize needs to be specified accordingly
 *
 * // Option B: High-Performance-Parallel Descriptor Calculation with SMILES string array
 * boolean result = Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
 *                      descriptors,
 *                      moleculeSmilesStrings,
 *                      matrix,
 *                      startIndex,
 *                      electronDonationModel,
 *                      isParallelCalculation,
 *                      nanPositionsListSynchronized
 * );
 * // you can also use Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(...) -> batchSize needs to be specified accordingly
 * }</pre>
 *
 * @author Manuel Schauer
 * @author Jonas Schaub
 * @author Achim Zielesny
 * @version 1.0.0
 */
public enum Descriptor {
    /*
     * Note for developers: for adding a new descriptor, go to "Add new descriptor information here!"
     * (also marked in the test class)
     */
    /**
     * Molecular weight, adds up the natural masses (weighted average of all known
     * isotopes of the particular element based on their natural abundances) of every atom
     * in the given molecule, it is NOT the exact mass.
     *
     * @see WeightDescriptor
     */
    MOLECULAR_WEIGHT(true, true, false, false, 1, "Molecular Weight"),
    /**
     * Wiener number, returns Wiener path number and Wiener polarity number.
     * Path number: sum of the distances between any two atoms in the molecule.<br>
     * Polarity number: number of pairs of atoms which are separated by exactly three bonds.<br>
     * Note: the CDK implementation counts all distances, not just those of carbon atoms or only carbon-carbon bonds.
     *
     * @see WienerNumbersDescriptor
     */
    WIENER_NUMBER(true, true, false, false, 2, "Wiener Number"),
    /**
     * Atom count, counts the number of all atoms in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT(true, true, false, false, 1, "Atom Count"),
    /**
     * Atom count heavy, counts the number of all heavy atoms in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_HEAVY(true, true, false, false, 1, "Atom Count"),
    /**
     * Atom count C, counts the number of all carbon atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_C(true, true, false, false, 1, "Atom Count C"),
    /**
     * Atom count H, counts the number of all hydrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_H(true, true, false, false, 1, "Atom Count H"),
    /**
     * Atom count N counts the number of all nitrogen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_N(true, true, false, false, 1, "Atom Count N"),
    /**
     * Atom count O, counts the number of all oxygen atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_O(true, true, false, false, 1, "Atom Count O"),
    /**
     * Atom count S, counts the number of all sulfur atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_S(true, true, false, false, 1, "Atom Count S"),
    /**
     * Atom count P, counts the number of all phosphorus atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_P(true, true, false, false, 1, "Atom Count P"),
    /**
     * Atom count F, counts the number of all fluorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_F(true, true, false, false, 1, "Atom Count F"),
    /**
     * Atom count Br, counts the number of all bromine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_BR(true, true, false, false, 1, "Atom Count Br"),
    /**
     * Atom count Cl, counts the number of all chlorine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_CL(true, true, false, false, 1, "Atom Count Cl"),
    /**
     * Atom count I, counts the number of all iodine atoms separately in the given molecule.
     *
     * @see AtomCountDescriptor
     */
    ATOM_COUNT_I(true, true, false, false, 1, "Atom Count I"),
    /**
     * Total bond count, counts the number of all bonds in a molecule, neglecting the order.
     * Double and triple bonds are counted as one bond.
     * Bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_ALL(true, true, false, false, 1, "Bond Count All"),
    /**
     * Bond count single, counts the number of single bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_SINGLE(true, true, false, false, 1, "Bond Count Single"),
    /**
     * Bond count double, counts the number of double bonds in a molecule.
     * No bonds to hydrogen atoms are counted.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_DOUBLE(true, true, false, false, 1, "Bond Count Double"),
    /**
     * Bond count triple, counts the number of triple bonds in a molecule.
     *
     * @see BondCountDescriptor
     */
    BOND_COUNT_TRIPLE(true, true, false, false, 1, "Bond Count Triple"),
    /**
     * H bond acceptor count, counts hydrogen bond acceptors based on a simplified PHACIR scheme.
     * It includes: Oxygen atoms with formal charge ≤ 0 (excluding: Aromatic ether oxygens and oxygens adjacent to nitrogen)
     * and nitrogen atoms with formal charge ≤ 0 (excluding: nitrogens adjacent to oxygen).
     *
     * @see HBondAcceptorCountDescriptor
     */
    H_BOND_ACCEPTOR_COUNT(true, true, false, false, 1, "H-Bond Acceptor Count"),
    /**
     * H bond donor count, counts hydrogen bond donors based on a simplified PHACIR classification.
     * It includes: OH groups where the oxygen has a formal charge ≥ 0 and NH groups where the nitrogen has a formal charge ≥ 0.
     *
     * @see HBondDonorCountDescriptor
     */
    H_BOND_DONOR_COUNT(true, true, false, false, 1, "H-Bond Donor Count"),
    /**
     * Aromatic atoms count, counts the number of aromatic atoms in a molecule.
     * Note: Requires that aromatic atoms in the molecule have already been detected and marked.
     *
     * @see AromaticAtomsCountDescriptor
     */
    AROMATIC_ATOMS_COUNT(true, true, false, false, 1, "Aromatic Atoms Count"),
    /**
     * Aromatic bonds count, counts the number of aromatic bonds in a molecule.
     * Note: Requires that aromatic bonds in the molecule have already been detected and marked.
     *
     * @see AromaticBondsCountDescriptor
     */
    AROMATIC_BONDS_COUNT(true, true, false, false, 1, "Aromatic Bonds Count"),
    /**
     * Rotatable bonds count, counts the number of rotatable bonds in a molecule.
     * A rotatable bond is defined as any single non-ring bond, where atoms on both sides
     * have at least two heavy-atom neighbors. Excluding terminal bonds.
     *
     * @see RotatableBondsCountDescriptor
     */
    ROTATABLE_BONDS_COUNT(true, true, false, false, 1, "Rotatable Bonds Count"),
    /**
     * Basic group count, returns the number of basic groups in a molecule.
     *
     * @see BasicGroupCountDescriptor
     */
    BASIC_GROUP_COUNT(true, true, false, false, 1, "Basic Group Count"),
    /**
     * Acidic group count, returns the number of acidic groups in a molecule.
     *
     * @see AcidicGroupCountDescriptor
     */
    ACIDIC_GROUP_COUNT(true, true, false, false, 1, "Acidic Group Count"),
    /**
     * TPSA descriptor, calculates the topological polar surface area (TPSA) of a molecule.
     * TPSA is the sum of the surface areas of polar atoms (typically oxygen and nitrogen)
     * and their attached hydrogens, based on a topological approximation (2D structure only).
     *
     * @see TPSADescriptor
     */
    TPSA(true, true, false, false, 1, "TPSA"),
    /**
     * Largest chain descriptor, calculates the number of atoms in the longest chain in the molecule.
     * This is a simple topological descriptor that provides a measure of molecular linearity.
     *
     * @see LargestChainDescriptor
     */
    LARGEST_CHAIN(true, true, false, false, 1, "Largest Chain"),
    /**
     * Longest aliphatic chain descriptor, calculates the number of atoms in the longest aliphatic chain.
     * This descriptor provides information about the maximum linear extent of non-aromatic
     * portions of the molecular structure, which relates to molecular shape properties.
     *
     * @see LongestAliphaticChainDescriptor
     */
    LONGEST_ALIPHATIC_CHAIN(true, true, false, false, 1, "Longest Aliphatic Chain"),
    /**
     * BCUT descriptor, calculates Burden matrix modified eigenvalues with different weighting schemes. Returns 6 values:<br>
     * 1. BCUTw-1l, BCUTw-2l ... - nlow lowest atom weighted BCUTS<br>
     * 2. BCUTw-1h, BCUTw-2h ... - nhigh highest atom weighted BCUTS<br>
     * 3. BCUTc-1l, BCUTc-2l ... - nlow lowest partial charge weighted BCUTS<br>
     * 4. BCUTc-1h, BCUTc-2h ... - nhigh highest partial charge weighted BCUTS<br>
     * 5. BCUTp-1l, BCUTp-2l ... - nlow lowest polarizability weighted BCUTS<br>
     * 6. BCUTp-1h, BCUTp-2h ... - nhigh highest polarizability weighted BCUTS<br>
     * Note: No array for one parameter is returned, just the highest and lowest numbers. (Default Parameters: nhigh = 1 and nlow = 1)
     *
     * @see BCUTDescriptor
     */
    BCUT(false, false, false, false, 6, "BCUT"),
    /**
     * Bond polarizability descriptor.
     * The BPolDescriptor calculates the bond polarizability of a molecule.
     * Bond polarizability is a simple sum of polarizability contributions from all bonds, based on bond types and involved atoms.
     * It provides a rough estimate of how easily the electron cloud in a molecule can be distorted,
     * which relates to intermolecular interactions, polarizability and refractive behavior.
     *
     * @see BPolDescriptor
     */
    B_POL(true, true, false, false, 1, "BPol"),
    /**
     * Rule of five descriptor, calculates the number of failures of Lipinski's Rule of Five.
     * The descriptor returns the number of violations (0-4).
     *
     * @see RuleOfFiveDescriptor
     */
    RULE_OF_FIVE(true, true, false, false, 1, "Rule of Five"),
    /**
     * FMF (Framework Match Fraction) descriptor, calculates the ratio of heavy atoms in
     * the framework to the total number of heavy atoms in the molecule.
     * This provides an indication of the proportion of the molecule that is part of the
     * scaffold or core structure, versus the proportion that is in side chains.
     *
     * @see FMFDescriptor
     */
    FMF(true, true, false, false, 1, "FMF"),
    /**
     * Fractional C SP3 descriptor, characterizes the non-flatness of a molecule by calculating
     * the fraction of sp3 hybridized carbon atoms over the total carbon count.
     * This provides information about the three-dimensionality and complexity
     * of a molecule, which relates to drug-likeness properties.
     *
     * @see FractionalCSP3Descriptor
     */
    FRACTIONAL_CSP3(true, true, false, false, 1, "Fractional CSP3"),
    /**
     * Hybridization ratio descriptor, calculates the ratio of sp3 carbons to sp2 carbons.
     * This provides valuable information about the three-dimensionality and flatness
     * of a molecule, which can be useful for predicting drug-like properties and
     * comparing structural characteristics.
     *
     * @see HybridizationRatioDescriptor
     */
    HYBRIDIZATION_RATIO(true, true, false, true, 1, "Hybridization Ratio"),
    /**
     * Kappa shape indices descriptor, calculates Kier and Hall kappa molecular shape indices.
     * These indices compare the molecular graph with minimal and maximal molecular graphs. Returns 3 values:<br>
     * 1. Kier1 - First kappa shape index<br>
     * 2. Kier2 - Second kappa shape index<br>
     * 3. Kier3 - Third kappa shape index<br>
     * Note: Hydrogens are ignored in the calculation.
     *
     * @see KappaShapeIndicesDescriptor
     */
    KAPPA_SHAPE_INDICES(false, true, false, true, 3, "Kappa Shape Indices"),
    /**
     * Petitjean number descriptor, calculates an index characterizing molecular graph topology.
     * This topological descriptor is based on the calculation of the graph eccentricity
     * and provides information about the molecular shape and branching pattern.
     *
     * @see PetitjeanNumberDescriptor
     */
    PETITJEAN_NUMBER(true, true, false, false, 1, "Petitjean Number"),
    /**
     * Spiro atom count descriptor, calculates the number of spiro atoms in a molecule.
     *
     * @see SpiroAtomCountDescriptor
     */
    SPIRO_ATOM_COUNT(true, true, false, false, 1, "Spiro Atom Count"),
    /**
     * VAdjMa descriptor, calculates the Vertex adjacency information (magnitude).
     * This is calculated as 1 + log2 m, where m is the number of heavy-heavy bonds.
     * If m is zero, then zero is returned.
     * This descriptor characterizes molecular complexity in terms of edge connectivity.
     *
     * @see VAdjMaDescriptor
     */
    V_ADJ_MAT(true, true, false, false, 1, "VAdjMa"),
    /**
     * Weighted path descriptor, evaluates the weighted path descriptors for a molecule.
     * Returns 5 values:<br>
     * 1. WTPT1 - molecular ID<br>
     * 2. WTPT2 - molecular ID / number of atoms<br>
     * 3. WTPT3 - sum of path lengths starting from heteroatoms<br>
     * 4. WTPT4 - sum of path lengths starting from oxygens<br>
     * 5. WTPT5 - sum of path lengths starting from nitrogens<br>
     * <p>
     * Note: This descriptor computes all paths which is an NP-hard problem, do not use it for complex molecules.
     * @see WeightedPathDescriptor
     */
    WEIGHTED_PATH(false, true, false, false, 5, "Weighted Path"),
    /**
     * Zagreb index descriptor, calculates the Zagreb index of a molecule.
     * The Zagreb index is the sum of the squares of atom degrees over all heavy atoms,
     * which provides information about the molecular complexity and topological structure.
     *
     * @see ZagrebIndexDescriptor
     */
    ZAGREB_INDEX(true, true, false, false, 1, "Zagreb Index"),
    /**
     * CarbonTypes descriptor, calculates the frequency of occurrence of 9 different types of carbon atoms. Returns 9 values:<br>
     * 1. C1SP1 - triply bound carbon bound to one other carbon<br>
     * 2. C2SP1 - triply bound carbon bound to two other carbons<br>
     * 3. C1SP2 - doubly bound carbon bound to one other carbon<br>
     * 4. C2SP2 - doubly bound carbon bound to two other carbons<br>
     * 5. C3SP2 - doubly bound carbon bound to three other carbons<br>
     * 6. C1SP3 - singly bound carbon bound to one other carbon<br>
     * 7. C2SP3 - singly bound carbon bound to two other carbons<br>
     * 8. C3SP3 - singly bound carbon bound to three other carbons<br>
     * 9. C4SP3 - singly bound carbon bound to four other carbons
     *
     * @see CarbonTypesDescriptor
     */
    CARBON_TYPES(true, true, false, false, 9, "Carbon Types"),
    /**
     * ALogP descriptor, calculates Ghose-Crippen LogP values, molar refractivity values
     * and ALogP squared values. Returns 3 values:<br>
     * 1. ALogP (logP value) is the Ghose-Crippen octanol-water partition coefficient.<br>
     * 2. ALogP² is the squared ALogP value.<br>
     * 3. Molar Refractivity (MR) measures the volume occupied by an atom or group of atoms.<br>
     *
     * @see ALOGPDescriptor
     */
    A_LOG_P(true, true, false, true, 3, "ALogP"),
    /**
     * XLogP descriptor, calculates logP based on the atom-type method called XLogP.
     * Requires all hydrogen's to be explicit.
     *
     * @see XLogPDescriptor
     */
    X_LOG_P(true, true, false, true, 1, "XLogP"),
    /**
     * JPlogP descriptor, calculates the octanol-water partition coefficient based on an atom contribution model.
     *
     * @see JPlogPDescriptor
     */
    JP_LOG_P(true, false, false, false, 1, "JPlogP"),
    /**
     * Mannhold LogP descriptor, calculates the octanol-water partition coefficient (logP) using the Mannhold method.
     * LogP describes the hydrophilicity or lipophilicity of a compound and is crucial for
     * predicting solubility, permeability, and bioavailability.
     *
     * @see MannholdLogPDescriptor
     */
    MANNHOLD_LOGP(true, true, false, false, 1, "Mannhold LogP"),
    /**
     * APol descriptor, calculates the sum of the atomic polarizabilities (including implicit hydrogens).
     *
     * @see APolDescriptor
     */
    A_POL(true, true, false, false, 1, "APol"),
    /**
     * Autocorrelation charge descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to charge distribution across the molecular structure.
     * This descriptor correlates atomic partial charges along the molecular topology
     * to characterize charge-related structural patterns in the molecule.
     * Returns 5 values representing charge autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorCharge
     */
    AUTOCORRELATION_CHARGE(true, false, false, false, 5, "Autocorrelation Charge"),
    /**
     * Autocorrelation mass descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to atomic mass distribution across the molecular structure.
     * This descriptor correlates atomic masses along the molecular topology
     * to characterize mass-related structural patterns in the molecule.
     * Returns 5 values representing mass autocorrelation at different topological distances.
     *
     * @see AutocorrelationDescriptorMass
     */
    AUTOCORRELATION_MASS(true, true, false, false, 5, "Autocorrelation Mass"),
    /**
     * Autocorrelation polarizability descriptor, calculates topological autocorrelation vectors
     * that capture patterns related to polarizability distribution across the molecular structure.
     * This descriptor correlates atomic polarizabilities along the molecular topology
     * to characterize polarizability-related structural patterns in the molecule.
     * Returns 5 values representing polarizability autocorrelation at different topological distances.
     * NOTE: Method is not validated in the CDK so not validated in this implementation as well
     *
     * @see AutocorrelationDescriptorPolarizability
     */
    AUTOCORRELATION_POLARIZABILITY(true, true, false, false, 5, "Autocorrelation Polarizability"),
    /**
     * Fragment complexity descriptor, calculates the complexity of a molecular system.
     * The complexity is defined as [Nilakantan, R. et al. Journal of chemical information and modeling. 2006. 46]:
     * C = abs(B^2 - A^2 + A) + H/100
     * where:
     * (C = complexity,
     * A = number of non-hydrogen atoms,
     * B = number of bonds,
     * H = number of heteroatoms,).
     * This provides a measure of structural complexity that correlates with synthetic accessibility.
     *
     * @see FragmentComplexityDescriptor
     */
    FRAGMENT_COMPLEXITY(true, true, false, false, 1, "Fragment Complexity"),
    /**
     * Chi chain descriptor, calculates the Kier + Hall chi chain indices of orders 3 through 7.
     * These values characterize a molecular graph based on its chain subgraphs.
     * Returns 10 values:<br>
     * 1. SCH-3 - Simple chain, order 3<br>
     * 2. SCH-4 - Simple chain, order 4<br>
     * 3. SCH-5 - Simple chain, order 5<br>
     * 4. SCH-6 - Simple chain, order 6<br>
     * 5. SCH-7 - Simple chain, order 7<br>
     * 6. VCH-3 - Valence chain, order 3<br>
     * 7. VCH-4 - Valence chain, order 4<br>
     * 8. VCH-5 - Valence chain, order 5<br>
     * 9. VCH-6 - Valence chain, order 6<br>
     * 10. VCH-7 - Valence chain, order 7
     *
     * @see ChiChainDescriptor
     */
    CHI_CHAIN(false, true, false, false, 10, "Chi Chain"),
    /**
     * Chi cluster descriptor, calculates Kier + Hall chi cluster indices of orders 3 through 6.
     * These values characterize a molecular graph based on its cluster subgraphs.
     * Returns 8 values:<br>
     * 1. SC-3 - Simple cluster, order 3<br>
     * 2. SC-4 - Simple cluster, order 4<br>
     * 3. SC-5 - Simple cluster, order 5<br>
     * 4. SC-6 - Simple cluster, order 6<br>
     * 5. VC-3 - Valence cluster, order 3<br>
     * 6. VC-4 - Valence cluster, order 4<br>
     * 7. VC-5 - Valence cluster, order 5<br>
     * 8. VC-6 - Valence cluster, order 6
     *
     * @see ChiClusterDescriptor
     */
    CHI_CLUSTER(false, true, false, false, 8, "Chi Cluster"),
    /**
     * Chi path cluster descriptor, calculates Kier + Hall chi path cluster indices of orders 4 through 6.
     * These values characterize a molecular graph based on its path cluster subgraphs.
     * Returns 6 values:<br>
     * 1. SPC-4 - Simple path cluster, order 4<br>
     * 2. SPC-5 - Simple path cluster, order 5<br>
     * 3. SPC-6 - Simple path cluster, order 6<br>
     * 4. VPC-4 - Valence path cluster, order 4<br>
     * 5. VPC-5 - Valence path cluster, order 5<br>
     * 6. VPC-6 - Valence path cluster, order 6
     *
     * @see ChiPathClusterDescriptor
     */
    CHI_PATH_CLUSTER(false, true, false, false, 6, "Chi Path Cluster"),
    /**
     * Chi path descriptor, calculates Kier + Hall chi path indices of orders 0 through 7.
     * These values characterize a molecular graph based on its path subgraphs.
     * Returns 16 values:<br>
     * 1.  SP-0 - Simple path, order 0<br>
     * 2.  SP-1 - Simple path, order 1<br>
     * 3.  SP-2 - Simple path, order 2<br>
     * 4.  SP-3 - Simple path, order 3<br>
     * 5.  SP-4 - Simple path, order 4<br>
     * 6.  SP-5 - Simple path, order 5<br>
     * 7.  SP-6 - Simple path, order 6<br>
     * 8.  SP-7 - Simple path, order 7<br>
     * 9.  VP-0 - Valence path, order 0<br>
     * 10. VP-1 - Valence path, order 1<br>
     * 11. VP-2 - Valence path, order 2<br>
     * 12. VP-3 - Valence path, order 3<br>
     * 13. VP-4 - Valence path, order 4<br>
     * 14. VP-5 - Valence path, order 5<br>
     * 15. VP-6 - Valence path, order 6<br>
     * 16. VP-7 - Valence path, order 7<br>
     *
     * @see ChiPathDescriptor
     */
    CHI_PATH(false, true, false, false, 16, "Chi Path"),
    /**
     * Fractional PSA descriptor, calculates the ratio of polar surface area to molecular weight.
     * This descriptor provides the polar surface area efficiency, which is the TPSADescriptor value divided by the
     * molecular weight, measured in square Angstroms per Dalton.
     *
     * @see FractionalPSADescriptor
     */
    FRACTIONAL_PSA(true, true, false, false, 1, "Fractional PSA"),
    /**
     * Largest pi system descriptor, calculates the number of atoms in the largest pi system.
     * This descriptor identifies the largest conjugated pi system within a molecule and
     * returns the count of atoms participating in it.
     *
     * @see LargestPiSystemDescriptor
     */
    LARGEST_PI_SYSTEM(true, true, false, false, 1, "Largest Pi System"),
    /**
     * Descriptor that calculates small ring information.
     * Returns 11 values:<br>
     * 1. nSmallRings - total number of small rings (of size 3 through 9)<br>
     * 2. nAromRings - total number of small aromatic rings<br>
     * 3. nRingBlocks - total number of distinct ring blocks<br>
     * 4. nAromBlocks - total number of aromatically connected components<br>
     * 5. nRings3 - total number of 3-membered rings<br>
     * 6. nRings4 - total number of 4-membered rings<br>
     * 7. nRings5 - total number of 5-membered rings<br>
     * 8. nRings6 - total number of 6-membered rings<br>
     * 9. nRings7 - total number of 7-membered rings<br>
     * 10. nRings8 - total number of 8-membered rings<br>
     * 11. nRings9 - total number of 9-membered rings<br>
     *
     * @see SmallRingDescriptor
     */
    SMALL_RING(true, true, false, false, 11, "Small Ring"),
    /**
     * Amino acid count descriptor, calculates the number of each amino acid in a molecule.
     * Returns 20 values, one for each of the 20 standard amino acids:
     * Alanine, Arginine, Asparagine, Aspartic acid, Cysteine, Glutamic acid, Glutamine,
     * Glycine, Histidine, Isoleucine, Leucine, Lysine, Methionine, Phenylalanine,
     * Proline, Serine, Threonine, Tryptophan, Tyrosine, and Valine.
     * This descriptor helps identify and quantify amino acid composition in peptides and proteins.
     *
     * @see AminoAcidCountDescriptor
     */
    AMINO_ACID_COUNT(false, true, false, false, 20, "Amino Acid Count"),
    /**
     * Kier-Hall SMARTS descriptor that calculates counts of functional groups and substructures
     * based on the Kier and Hall SMARTS patterns, used for QSAR modeling and molecular characterization.
     * Note: This descriptor provides 79 values representing different molecular fragments.
     *
     * @see KierHallSmartsDescriptor
     */
    KIER_HALL_SMARTS(true, true, false, false, 79, "Kier Hall SMARTS"),
    /**
     * Eccentric connectivity index descriptor, calculates a topological descriptor that combines
     * distance and adjacency information.
     * It is defined as the sum of the products of eccentricity and vertex degree for each atom.
     * This index provides information about the distribution of atoms in the molecular structure
     * and helps characterize molecular complexity, branching, and overall shape.
     *
     * @see EccentricConnectivityIndexDescriptor
     */
    ECCENTRIC_CONNECTIVITY_INDEX(true, true, false, false, 1, "Eccentric Connectivity Index"),
    /**
     * MDE descriptor, calculates molecular distance edge descriptors for carbon, oxygen and nitrogen atoms.
     * These descriptors encode information about the connectivity and distance of atoms of specific types
     * and hybridization states in the molecular graph. Returns 19 values representing various molecular
     * distance edge counts between different types of carbon, oxygen, and nitrogen atoms:<br>
     * MDEC-11<br>
     * MDEC-12<br>
     * MDEC-13<br>
     * MDEC-14<br>
     * MDEC-22<br>
     * MDEC-23<br>
     * MDEC-24<br>
     * MDEC-33<br>
     * MDEC-34<br>
     * MDEC-44<br>
     * MDEO-11<br>
     * MDEO-12<br>
     * MDEO-22<br>
     * MDEN-11<br>
     * MDEN-12<br>
     * MDEN-13<br>
     * MDEN-22<br>
     * MDEN-23<br>
     * MDEN-33<br>
     *
     * @see MDEDescriptor
     */
    MDE(true, true, false, false, 19, "MDE"),
    /**
     * VABC descriptor, calculates the volume descriptor using the van der Waals volume calculation approach.
     * This descriptor estimates molecular volume based on atom contributions, considering bond types
     * and atomic properties, providing insights into molecular size and steric properties.
     *
     * @see VABCDescriptor
     */
    VABC(true, true, false, false, 1, "VABC"),
    /**
     * PubChem fingerprinter, generates a 881-bit binary fingerprint based on PubChem's substructure keys.
     * This fingerprint encodes the presence or absence of specific substructural features
     * and is useful for similarity searching and chemical space analysis.
     *
     * @see PubchemFingerprinter
     */
    PUBCHEM_FINGERPRINTER(false, true, true, false, 881, "PubChem Fingerprinter"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_0(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter ECFP-0"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 0.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_0(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter FCFP-0"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_2(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter ECFP-2"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 2.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_2(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter FCFP-2"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_4(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter ECFP-4"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 4.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_4(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter FCFP-4"),
    /**
     * Circular fingerprinter, generates an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_ECFP_6(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter ECFP-6"),
    /**
     * Circular fingerprinter, generates a functional class version of an extended-connectivity fingerprint with a path diameter of 6.
     *
     * @see CircularFingerprinter
     */
    CIRCULAR_FINGERPRINTER_FCFP_6(true, true, true, false, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, "Circular Fingerprinter FCFP-6"),
    /**
     * MACCS fingerprinter, generates a 166-bit binary fingerprint based on the MACCS structural keys.
     *
     * @see MACCSFingerprinter
     */
    MACCS_FINGERPRINTER(true, true, true, false, 166, "MACCS Fingerprinter");

    // Add new descriptor information here!

    /**
     * Indicates whether this descriptor is quickly calculable.
     * Fast descriptors have lower computational complexity and can be calculated efficiently,
     * while slow descriptors require more intensive calculations (e.g., CHI descriptors).
     */
    private final boolean isFast;

    /**
     * Indicates whether this descriptor is safe to calculate.
     * Safe descriptors are reliable and validated, while unsafe descriptors may produce
     * inconsistent results or NaN values under certain conditions.
     */
    private final boolean isSafe;

    /**
     * Indicates whether the enum entry is a fingerprint or not.
     */
    private final boolean isFingerprint;

    /**
     * Indicates whether this descriptor needs explicit hydrogens.
     */
    private final boolean needsExplicitHydrogens;

    /**
     * The number of components calculated by this descriptor.
     */
    private int descriptorComponentNumber;

    /**
     * The human-readable name of this descriptor for output purposes.
     */
    private final String name;

    /**
     * Constructs a Descriptor with the given field values.
     *
     * @param isFast true if the descriptor is quickly calculable, false if it requires
     *               more intensive computation
     * @param isSafe true if the descriptor is safe and reliable, false if it may produce
     *               inconsistent results or NaN values
     * @param isFingerprint true if the "descriptor" is actually a fingerprint, false otherwise
     * @param needsExplicitHydrogens true if the descriptor requires explicit hydrogens, false otherwise
     * @param descriptorComponentNumber the number of components calculated by this descriptor
     * @param name the human-readable name of this descriptor for output purposes
     */
    Descriptor(boolean isFast, boolean isSafe, boolean isFingerprint, boolean needsExplicitHydrogens, int descriptorComponentNumber, String name) {
        this.isFast = isFast;
        this.isSafe = isSafe;
        this.isFingerprint = isFingerprint;
        this.needsExplicitHydrogens = needsExplicitHydrogens;
        this.descriptorComponentNumber = descriptorComponentNumber;
        this.name = name;
    }

    /**
     * Returns whether this descriptor can be calculated quickly. This classification is based on a performance benchmark
     * of 10,000 molecules: descriptors computed in under 1 second for 10,000 molecules are considered fast.
     *
     * @return true if the descriptor has low computational complexity and can be
     *         calculated efficiently, false if it requires intensive computation
     */
    public boolean isFast() {
        return this.isFast;
    }

    /**
     * Returns whether this descriptor is safe to calculate.
     *
     * @return true if the descriptor is reliable and validated, false if it may
     *         produce inconsistent results or NaN values under certain conditions
     */
    public boolean isSafe() {
        return this.isSafe;
    }

    /**
     * Returns whether this descriptor calculates a fingerprint.
     *
     * @return true if the enum entry calculates a fingerprint, false otherwise
     */
    public boolean isFingerprint() {
        return this.isFingerprint;
    }

    /**
     * Returns whether this descriptor needs explicit hydrogens.
     *
     * @return true if the descriptor requires explicit hydrogens, false otherwise
     */
    public boolean needsExplicitHydrogens() {
        return this.needsExplicitHydrogens;
    }

    /**
     * Returns the number of components calculated by this descriptor.
     *
     * @return the number of components calculated by this descriptor
     */
    public int getDescriptorComponentNumber() {
        return this.descriptorComponentNumber;
    }

    /**
     * Returns the human-readable name of this descriptor.
     *
     * @return the human-readable name of this descriptor
     */
    public String getName() {
        return this.name;
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
     * @param atomContainer the molecule to calculate descriptors for (IS NOT CHANGED);
     *                        must have aromaticity already perceived if required by the descriptor
     * @param vector         the result vector to store calculated values (MAY BE CHANGED);
     *                        values are stored starting at startIndex
     * @param startIndex     the starting index in vector where results should be written;
     *                        subsequent values are written to startIndex + 1, startIndex + 2, etc.
     * @return Returns true if the calculation was successful and results were stored in the vector, false otherwise.
     */
    private boolean calculate(IAtomContainer atomContainer, float[] vector, int startIndex) throws InterruptedException {
        boolean success = false;
        try {
            if (this.isFingerprint()) {
                return this.calculateFingerprintFromPool(atomContainer, vector, startIndex);
            } else {
                IMolecularDescriptor cdkDescriptor = Descriptor.descriptorToCdkObjectMap.get(this);
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
                        for (int i = 0; i < this.descriptorComponentNumber; i++) {
                            vector[startIndex + i] = (float) arrayResult.get(i);
                        }
                        success = true;
                    } else if (result instanceof IntegerArrayResult) {
                        IntegerArrayResult arrayResult = (IntegerArrayResult) result;
                        for (int i = 0; i < this.descriptorComponentNumber; i++) {
                            vector[startIndex + i] = (float) arrayResult.get(i);
                        }
                        success = true;
                    } else {
                        LOGGER.log(Level.WARNING, () -> "Unexpected result type " + result.getClass().getName() + " for descriptor " + this.name() + ". Expected DoubleResult, IntegerResult, DoubleArrayResult or IntegerArrayResult.");
                    }
                } else {
                    LOGGER.log(Level.WARNING, () -> "Descriptor " + this.name() + " not found in descriptorToCdkObjectMap. This should not happen.");
                }
            }
        } catch (InterruptedException exception) {
            LOGGER.log(Level.WARNING, exception, () -> "Interrupted while waiting for fingerprinter instance from pool. This should not happen.");
            throw exception;
        } catch (Exception exception) {
            LOGGER.log(Level.SEVERE, exception, () -> "Descriptor " + this.name() + " failed to calculate.");
            return false;
        }
        return success;
    }

    /**
     * Logger of this class.
     * TODO for CDK integration: must be replaced with CDK ILoggingTool instance.
     */
    static final Logger LOGGER = Logger.getLogger(Descriptor.class.getName());

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
     * Default size used for all circular fingerprint "descriptors".
     * Note: This constant must remain static final because it is required during
     * the static initialization of the enum constants above. Enum constants are
     * instantiated before static variables are initialized, so using a non-final
     * variable would result in a default value of 0.
     */
    private static final int CIRCULAR_FINGERPRINT_DEFAULT_SIZE = 1024;

    /**
     * Size used for all circular fingerprint "descriptors". Default value is 1024.
     * Can be changed via the synchronized {@link #setCircularFingerprintSize(int)} method.
     */
    private static volatile int circularFingerprintSize = CIRCULAR_FINGERPRINT_DEFAULT_SIZE;
    /*
     * Static initializer block to populate the descriptorToCdkObjectMap and initialize the fingerprint pool.
     * Note: We use the map and initialize it here (instead of giving each descriptor constant an instance field)
     * to be able to do error handling, and because IFingerprinter and IMolecularDescriptor are different object types.
     */
    static {
        try {
            // MOLECULAR_WEIGHT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.MOLECULAR_WEIGHT, new WeightDescriptor());

            // WIENER_NUMBER
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.WIENER_NUMBER, new WienerNumbersDescriptor());

            // ATOM_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT, new AtomCountDescriptor());

            // ATOM_COUNT_HEAVY
            AtomCountDescriptor atomCountHeavyDescriptor = new AtomCountDescriptor();
            atomCountHeavyDescriptor.setParameters(new Object[] {"#"}); // set heavy atom count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_HEAVY, atomCountHeavyDescriptor);

            // ATOM_COUNT_C
            AtomCountDescriptor atomCountCDescriptor = new AtomCountDescriptor();
            atomCountCDescriptor.setParameters(new Object[] {"C"}); // set carbon as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_C, atomCountCDescriptor);

            // ATOM_COUNT_H
            AtomCountDescriptor atomCountHDescriptor = new AtomCountDescriptor();
            atomCountHDescriptor.setParameters(new Object[] {"H"}); // set hydrogen as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_H, atomCountHDescriptor);

            // ATOM_COUNT_N
            AtomCountDescriptor atomCountNDescriptor = new AtomCountDescriptor();
            atomCountNDescriptor.setParameters(new Object[] {"N"}); // set nitrogen as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_N, atomCountNDescriptor);

            // ATOM_COUNT_O
            AtomCountDescriptor atomCountODescriptor = new AtomCountDescriptor();
            atomCountODescriptor.setParameters(new Object[] {"O"}); // set oxygen as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_O, atomCountODescriptor);

            // ATOM_COUNT_S
            AtomCountDescriptor atomCountSDescriptor = new AtomCountDescriptor();
            atomCountSDescriptor.setParameters(new Object[] {"S"}); // set sulfur as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_S, atomCountSDescriptor);

            // ATOM_COUNT_P
            AtomCountDescriptor atomCountPDescriptor = new AtomCountDescriptor();
            atomCountPDescriptor.setParameters(new Object[] {"P"}); // set phosphorus as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_P, atomCountPDescriptor);

            // ATOM_COUNT_F
            AtomCountDescriptor atomCountFDescriptor = new AtomCountDescriptor();
            atomCountFDescriptor.setParameters(new Object[] {"F"}); // set fluorine as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_F, atomCountFDescriptor);

            // ATOM_COUNT_BR
            AtomCountDescriptor atomCountBrDescriptor = new AtomCountDescriptor();
            atomCountBrDescriptor.setParameters(new Object[] {"Br"}); // set bromine as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_BR, atomCountBrDescriptor);

            // ATOM_COUNT_CL
            AtomCountDescriptor atomCountClDescriptor = new AtomCountDescriptor();
            atomCountClDescriptor.setParameters(new Object[] {"Cl"}); // set chlorine as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_CL, atomCountClDescriptor);

            // ATOM_COUNT_I
            AtomCountDescriptor atomCountIDescriptor = new AtomCountDescriptor();
            atomCountIDescriptor.setParameters(new Object[] {"I"}); // set iodine as the atom to count
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ATOM_COUNT_I, atomCountIDescriptor);

            // H_BOND_ACCEPTOR_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.H_BOND_ACCEPTOR_COUNT, new HBondAcceptorCountDescriptor());

            // H_BOND_DONOR_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.H_BOND_DONOR_COUNT, new HBondDonorCountDescriptor());

            // TPSA
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.TPSA, new TPSADescriptor());

            // LARGEST_CHAIN
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.LARGEST_CHAIN, new LargestChainDescriptor());

            // LONGEST_ALIPHATIC_CHAIN
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.LONGEST_ALIPHATIC_CHAIN, new LongestAliphaticChainDescriptor());

            // MANNHOLD_LOGP
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.MANNHOLD_LOGP, new MannholdLogPDescriptor());

            // BCUT
            BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
            bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BCUT, bcutDescriptor);

            // BOND_COUNT_ALL
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_ALL, new BondCountDescriptor());

            // BOND_COUNT_SINGLE
            BondCountDescriptor bondCountSingleDescriptor = new BondCountDescriptor();
            bondCountSingleDescriptor.setParameters(new Object[]{"s"});
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_SINGLE, bondCountSingleDescriptor);

            // BOND_COUNT_DOUBLE
            BondCountDescriptor bondCountDoubleDescriptor = new BondCountDescriptor();
            bondCountDoubleDescriptor.setParameters(new Object[]{"d"});
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_DOUBLE, bondCountDoubleDescriptor);

            // BOND_COUNT_TRIPLE
            BondCountDescriptor bondCountTripleDescriptor = new BondCountDescriptor();
            bondCountTripleDescriptor.setParameters(new Object[]{"t"});
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BOND_COUNT_TRIPLE, bondCountTripleDescriptor);

            // B_POL
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.B_POL, new BPolDescriptor());

            // RULE_OF_FIVE
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.RULE_OF_FIVE, new RuleOfFiveDescriptor());

            // AROMATIC_ATOMS_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AROMATIC_ATOMS_COUNT, new AromaticAtomsCountDescriptor());

            // AROMATIC_BONDS_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AROMATIC_BONDS_COUNT, new AromaticBondsCountDescriptor());

            // ROTATABLE_BONDS_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ROTATABLE_BONDS_COUNT, new RotatableBondsCountDescriptor());

            // FMF
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.FMF, new FMFDescriptor());

            // FRACTIONAL_CSP3
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.FRACTIONAL_CSP3, new FractionalCSP3Descriptor());

            // HYBRIDIZATION_RATIO
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.HYBRIDIZATION_RATIO, new HybridizationRatioDescriptor());

            // KAPPA_SHAPE_INDICES
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.KAPPA_SHAPE_INDICES, new KappaShapeIndicesDescriptor());

            // PETITJEAN_NUMBER
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.PETITJEAN_NUMBER, new PetitjeanNumberDescriptor());

            // SPIRO_ATOM_COUNT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.SPIRO_ATOM_COUNT, new SpiroAtomCountDescriptor());

            // V_ADJ_MAT
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.V_ADJ_MAT, new VAdjMaDescriptor());

            // WEIGHTED_PATH
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.WEIGHTED_PATH, new WeightedPathDescriptor());

            // ZAGREB_INDEX
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ZAGREB_INDEX, new ZagrebIndexDescriptor());

            // CARBON_TYPES
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.CARBON_TYPES, new CarbonTypesDescriptor());

            // A_LOG_P
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.A_LOG_P, new ALOGPDescriptor());

            // X_LOG_P
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.X_LOG_P, new XLogPDescriptor());

            // JP_LOG_P
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.JP_LOG_P, new JPlogPDescriptor());

            // A_POL
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.A_POL, new APolDescriptor());

            // AUTOCORRELATION_CHARGE
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_CHARGE, new AutocorrelationDescriptorCharge());

            // AUTOCORRELATION_MASS
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_MASS, new AutocorrelationDescriptorMass());

            // AUTOCORRELATION_POLARIZABILITY
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AUTOCORRELATION_POLARIZABILITY, new AutocorrelationDescriptorPolarizability());

            // FRAGMENT_COMPLEXITY
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.FRAGMENT_COMPLEXITY, new FragmentComplexityDescriptor());

            // CHI_CHAIN
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.CHI_CHAIN, new ChiChainDescriptor());

            // CHI_CLUSTER
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.CHI_CLUSTER, new ChiClusterDescriptor());

            // CHI_PATH_CLUSTER
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.CHI_PATH_CLUSTER, new ChiPathClusterDescriptor());

            // CHI_PATH
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.CHI_PATH, new ChiPathDescriptor());

            // FRACTIONAL_PSA
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.FRACTIONAL_PSA, new FractionalPSADescriptor());

            // LARGEST_PI_SYSTEM
            LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
            largestPiSystemDescriptor.setParameters(new Object[] {false}); //do not check aromaticity again
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.LARGEST_PI_SYSTEM, largestPiSystemDescriptor);

            // SMALL_RING
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.SMALL_RING, new SmallRingDescriptor());

            // Basic Group Count
            BasicGroupCountDescriptor basicGroupCountDescriptor = new BasicGroupCountDescriptor();
            basicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.BASIC_GROUP_COUNT, basicGroupCountDescriptor);

            // Acidic Group Count
            AcidicGroupCountDescriptor acidicGroupCountDescriptor = new AcidicGroupCountDescriptor();
            acidicGroupCountDescriptor.initialise(SilentChemObjectBuilder.getInstance());
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ACIDIC_GROUP_COUNT, acidicGroupCountDescriptor);

            // AminoAcidCount
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.AMINO_ACID_COUNT, new AminoAcidCountDescriptor());

            // Kier-Hall SMARTS
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.KIER_HALL_SMARTS, new KierHallSmartsDescriptor());

            // ECCENTRIC_CONNECTIVITY_INDEX
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.ECCENTRIC_CONNECTIVITY_INDEX, new EccentricConnectivityIndexDescriptor());

            // MDE
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.MDE, new MDEDescriptor());

            // VABC
            Descriptor.descriptorToCdkObjectMap.put(Descriptor.VABC, new VABCDescriptor());

            // Add new descriptor information here!

            // Initialize fingerprint pools with configurable pool size (default: 4)
            Descriptor.initializeFingerprintPools();

        } catch (Exception exception) {
            throw new UnsupportedOperationException("Failed to initialize descriptors, this should never happen. ", exception);
        }
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
        Descriptor.fingerprintPoolMap.clear();

        // PUBCHEM_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> pubchemPool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!pubchemPool.offer(new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()))){
                LOGGER.log(Level.WARNING, () -> "Failed to add PubchemFingerprinter instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.PUBCHEM_FINGERPRINTER, pubchemPool);

        // CIRCULAR_FINGERPRINTER_ECFP_0 Pool
        BlockingQueue<IFingerprinter> ecfp0Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!ecfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP0, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter ECFP0 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0, ecfp0Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_0 Pool
        BlockingQueue<IFingerprinter> fcfp0Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!fcfp0Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP0, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter FCFP0 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0, fcfp0Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_2 Pool
        BlockingQueue<IFingerprinter> ecfp2Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!ecfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP2, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter ECFP2 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2, ecfp2Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_2 Pool
        BlockingQueue<IFingerprinter> fcfp2Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!fcfp2Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP2, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter FCFP2 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2, fcfp2Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_4 Pool
        BlockingQueue<IFingerprinter> ecfp4Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!ecfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter ECFP4 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, ecfp4Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_4 Pool
        BlockingQueue<IFingerprinter> fcfp4Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!fcfp4Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP4, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter FCFP4 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4, fcfp4Pool);

        // CIRCULAR_FINGERPRINTER_ECFP_6 Pool
        BlockingQueue<IFingerprinter> ecfp6Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!ecfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter ECFP6 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6, ecfp6Pool);

        // CIRCULAR_FINGERPRINTER_FCFP_6 Pool
        BlockingQueue<IFingerprinter> fcfp6Pool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!fcfp6Pool.offer(new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6, Descriptor.circularFingerprintSize))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add CircularFingerprinter FCFP6 instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6, fcfp6Pool);

        // MACCS_FINGERPRINTER Pool
        BlockingQueue<IFingerprinter> maccsPool = new LinkedBlockingQueue<>(Descriptor.fingerprintPoolSize);
        for (int i = 0; i < Descriptor.fingerprintPoolSize; i++) {
            if (!maccsPool.offer(new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()))) {
                LOGGER.log(Level.WARNING, () -> "Failed to add MACCSFingerprinter instance to pool. This should not happen.");
                success = false;
            }
        }
        Descriptor.fingerprintPoolMap.put(Descriptor.MACCS_FINGERPRINTER, maccsPool);

        return success;
    }

    /**
     * Sets fingerprint component values in vector.
     *
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
    private boolean calculateFingerprintFromPool(IAtomContainer atomContainer, float[] vector, int startIndex)
            throws InterruptedException {
        IFingerprinter fingerprinter = null;
        boolean success;
        try {
            BlockingQueue<IFingerprinter> pool = Descriptor.fingerprintPoolMap.get(this);
            // take() blocks until a fingerprinter is available in the pool; throws InterruptedException
            // if the thread is interrupted while waiting (or already interrupted). InterruptedException is
            // neither CDKException nor RuntimeException, so it is not caught below and propagates naturally.
            fingerprinter = pool.take();

            IBitFingerprint bitFingerprint = fingerprinter.getBitFingerprint(atomContainer);
            for (int i = 0; i < this.getDescriptorComponentNumber(); i++) {
                vector[startIndex + i] = bitFingerprint.get(i) ? 1.0f : 0.0f;
            }
            success = true;
        } catch (CDKException | RuntimeException exception) {
            LOGGER.log(Level.WARNING, exception, () -> "Failed to calculate: " + this.name());
            return false;
        } finally {
            // Only return the fingerprinter if it was successfully taken from the pool.
            // If take() threw InterruptedException, fingerprinter is still null and offer(null)
            // would throw NullPointerException.
            if (fingerprinter != null && !Descriptor.fingerprintPoolMap.get(this).offer(fingerprinter)){
                LOGGER.log(Level.WARNING, () -> "Failed to return fingerprinter to pool: " + this.name());
                success = false;
            }

        }
        return success;
    }

    /**
     * Returns all available descriptors.
     *
     * @return All available descriptors
     */
    public static Descriptor[] getAllDescriptors(){
        return Descriptor.values();
    }

    /**
     * Returns all available fingerprint descriptors.
     *
     * @return All available fingerprint descriptors
     */
    public static Descriptor[] getAllFingerprints() {
        return Arrays.stream(Descriptor.values())
                .filter(Descriptor::isFingerprint)
                .toArray(Descriptor[]::new);
    }

    /**
     * Returns specified available descriptors.
     * Note: Fast, safe, and non-fingerprint descriptors are automatically included by default.
     *
     * @param isSlowlyCalculableDescriptorInclusion True: Slowly calculable descriptors are included in the result, false: Otherwise.
     * @param isUnsafeDescriptorInclusion True: Unsafe descriptors are included in the result, false: Unsafe descriptors
     *                                     are excluded from the result, this does not mean no NaN's can be produced.
     * @param isFingerprintAsDescriptorInclusion True: Fingerprint descriptors are included in the result, false: Fingerprint descriptors are excluded.
     * @return Specified descriptors
     */
    public static Descriptor[] getSpecifiedDescriptors(
            boolean isSlowlyCalculableDescriptorInclusion,
            boolean isUnsafeDescriptorInclusion,
            boolean isFingerprintAsDescriptorInclusion
    ) {
        // Initialize ArrayList with maximum possible capacity to avoid internal resizing during element addition
        List<Descriptor> result = new ArrayList<>(Descriptor.values().length);

        for (Descriptor descriptor : Descriptor.values()) {
            boolean includeDescriptor = true;

            // Exclude slow descriptors if not requested
            if (!descriptor.isFast() && !isSlowlyCalculableDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude unsafe descriptors if not requested
            if (!descriptor.isSafe() && !isUnsafeDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Exclude fingerprint descriptors if not requested
            if (descriptor.isFingerprint() && !isFingerprintAsDescriptorInclusion) {
                includeDescriptor = false;
            }

            // Only fast, safe, and non-fingerprint descriptors remain if all flags are false
            if (includeDescriptor) {
                result.add(descriptor);
            }
        }

        return result.toArray(new Descriptor[0]);
    }

    /**
     * Returns sum of number of calculated components of an array of defined descriptors.
     *
     * @param descriptors Array of descriptors (IS NOT CHANGED); if it is empty, 0 is returned;
     *                     may not be null (the array or one of its elements)
     * @return Sum of number of calculated components of array of descriptors
     */
    public static int getNumberOfComponents(
        Descriptor[] descriptors
    ) {
        // Checks
        if (descriptors == null) {
            throw new NullPointerException("Descriptor.getNumberOfComponents: descriptor is null.");
        }
        if (descriptors.length == 0) {
            return 0;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor == null) {
                throw new NullPointerException("Descriptor.getNumberOfComponents: At least one descriptor in descriptors is null.");
            }
        }

        int totalNumberOfComponents = 0;
        for (Descriptor descriptor : descriptors) {
            totalNumberOfComponents += descriptor.getDescriptorComponentNumber();
        }
        return totalNumberOfComponents;
    }


    /**
     * Returns the descriptor and component index for a given position in a descriptor array.
     * <p>
     * Example: Given descriptors [MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)]:
     * <ul>
     *   <li>Index 0 → MOLECULAR_WEIGHT, component 0</li>
     *   <li>Index 1 → BCUT, component 0</li>
     *   <li>...</li>
     *   <li>Index 6 → BCUT, component 5</li>
     *   <li>Index 7 → WIENER_NUMBER, component 0</li>
     *   <li>Index 8 → WIENER_NUMBER, component 1</li>
     * </ul>
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A two-element int array: [descriptor index in descriptors, component index within that descriptor]
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static int[] getDescriptorAndComponentIndex(Descriptor[] descriptors, int anIndex)
            throws IllegalArgumentException {
        // Checks
        if (descriptors == null) {
            throw new NullPointerException(
                    "Descriptor.getDescriptorAndComponentIndex: descriptors must not be null."
            );
        }
        if (descriptors.length == 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: descriptors must not be null or empty."
            );
        }
        if (anIndex < 0) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: anIndex must be >= 0."
            );
        }
        if (anIndex >= Descriptor.getNumberOfComponents(descriptors)) {
            throw new IllegalArgumentException(
                    "Descriptor.getDescriptorAndComponentIndex: anIndex (" + anIndex +
                            ") exceeds total component count (" + Descriptor.getNumberOfComponents(descriptors) + ")."
            );
        }

        // Calculate cumulative component counts and find the target descriptor
        int cumulativeCount = 0;
        int[] result = new int[2];
        for (int i = 0; i < descriptors.length; i++) {
            if (descriptors[i] == null) {
                throw new NullPointerException(
                        "Descriptor.getDescriptorAndComponentIndex: descriptors contains null element at index " + i
                );
            }

            int componentCount = descriptors[i].getDescriptorComponentNumber();
            int nextCumulativeCount = cumulativeCount + componentCount;

            // Check if anIndex falls within this descriptor's range
            if (anIndex < nextCumulativeCount) {
                result[0] = i;
                result[1] = anIndex - cumulativeCount;
                break;
            }

            cumulativeCount = nextCumulativeCount;
        }
        return result;
    }

    /**
     * Returns the descriptor and its component name for a given position in a descriptor array.
     * <p>
     * This is a convenience method that extends {@link #getDescriptorAndComponentIndex} by also
     * providing human-readable names for the descriptor and component.
     *
     * @param descriptors Array of descriptors (must not be null or empty; its elements must not be null)
     * @param anIndex The index in the descriptor array (must be &gt;= 0 and &lt; total component count)
     * @return A String array: [descriptor name, component index as string]
     * @throws IllegalArgumentException if descriptors is empty, anIndex is negative,
     *                                  or anIndex exceeds the total number of components
     */
    public static String[] getDescriptorAndComponentInfo(Descriptor[] descriptors, int anIndex)
            throws IllegalArgumentException {
        int[] indices = Descriptor.getDescriptorAndComponentIndex(descriptors, anIndex);
        Descriptor descriptor = descriptors[indices[0]];

        return new String[] {
                descriptor.getName(),
                String.valueOf(indices[1])
        };
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
            Descriptor.fingerprintPoolSize = aPoolSize;
            Descriptor.initializeFingerprintPools();
            return true;
        } catch (Exception exception){
            LOGGER.log(Level.SEVERE,exception , () ->"Failed to set fingerprint pool size and reinitialize pools: " + exception.getMessage());
            return false;
        }

    }

    /**
     * Returns the current pool size for fingerprinter instances.
     *
     * @return Current fingerprint pool size
     */
    public static synchronized int getFingerprintPoolSize() {
        return Descriptor.fingerprintPoolSize;
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
            Descriptor.circularFingerprintSize = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6.descriptorComponentNumber = aSize;
            Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6.descriptorComponentNumber = aSize;
            Descriptor.initializeFingerprintPools();
            return true;
        } catch (Exception exception) {
            LOGGER.log(Level.SEVERE,exception , () -> "Failed to set circular fingerprint size and reinitialize pools: " + exception.getMessage());
            return false;
        }
    }

    /**
     * Returns the current size used for circular fingerprints.
     *
     * @return Current circular fingerprint size
     */
    public static synchronized int getCircularFingerprintSize() {
        return Descriptor.circularFingerprintSize;
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
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor, IAtomContainer molecule)  {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }

        float[] result = new float[descriptor.getDescriptorComponentNumber()];
        try {
            boolean success;
            if (descriptor.needsExplicitHydrogens()) {
                IAtomContainer moleculeWithExplicitHydrogens = Descriptor.createMoleculeWithExplicitHydrogens(molecule);
                success = descriptor.calculate(moleculeWithExplicitHydrogens, result, 0);
            } else {
                success = descriptor.calculate(molecule, result, 0);
            }
            if (!success) {
                Arrays.fill(result, Float.NaN);
            }
        } catch (InterruptedException exception) {
            Thread.currentThread().interrupt();
            Arrays.fill(result, Float.NaN);
            Descriptor.LOGGER.log(Level.WARNING,exception , () -> "Interrupted while calculating descriptor " + descriptor.getName());
        } catch (Exception exception) {
            // Fill with NaN on failure and log the error
            Arrays.fill(result, Float.NaN);
            Descriptor.LOGGER.log(Level.WARNING, String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
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
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor,
                                              IAtomContainer molecule,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: descriptor must not be null.");
        }
        if (molecule == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: molecule must not be null.");
        }
        if (molecule.isEmpty()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: electronDonationModel must not be null.");
        }

        try {
            Descriptor.setAromaticity(molecule, electronDonationModel);
            return Descriptor.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            Descriptor.LOGGER.log(Level.WARNING, String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
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
     *         {@link #getDescriptorComponentNumber()}. Contains NaN values if calculation fails.
     */
    public static float[] calculateDescriptor(Descriptor descriptor,
                                              String smilesString,
                                              ElectronDonation electronDonationModel
    )   {
        // Checks
        if (descriptor == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: descriptor must not be null.");
        }
        if (smilesString == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: smilesString must not be null.");
        }
        if (smilesString.isBlank()) {
            return new float[0];
        }
        if (electronDonationModel == null) {
            throw new NullPointerException("Descriptor.calculateDescriptor: electronDonationModel must not be null.");
        }

        try {
            IAtomContainer molecule = Descriptor.SMILES_PARSER.parseSmiles(smilesString);
            Descriptor.setAromaticity(molecule, electronDonationModel);
            return Descriptor.calculateDescriptor(descriptor, molecule);
        } catch (Exception exception) {
            float[] result = new float[descriptor.getDescriptorComponentNumber()];
            Arrays.fill(result, Float.NaN);
            Descriptor.LOGGER.log(Level.WARNING, String.format("Failed to calculate descriptor %s: %s", descriptor.getName(), exception.getMessage()), exception);
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
        final String methodName = "setDescriptorsForMoleculesByMoleculeBatchParallelization";
        if (!Descriptor.validateDescriptors(descriptors, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given descriptor array is empty, calculation aborted.", methodName);
            return true;
        }
        if (!Descriptor.validateAtomContainerArray(atomContainerArray, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given atom container array is empty, calculation aborted.", methodName);
            return true;
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        Descriptor.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

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
                            boolean success = Descriptor.setDescriptorsForSingleMolecule(
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
                            Descriptor.LOGGER.log(
                                    Level.WARNING,
                                    String.format("Descriptor.setDescriptorsForMoleculesByBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    //TODO: why is there no try-catch for every single descriptor calculation necessary here? Same in the methods below
                    if (!Descriptor.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, nanPositionsList)) {
                        hasNaN.set(true);
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            Descriptor.LOGGER.log(
                    Level.WARNING,
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
        if (!Descriptor.validateDescriptors(descriptors, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given descriptor array is empty, calculation aborted.", methodName);
            return true;
        }
        if (!Descriptor.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given SMILES string array is empty, calculation aborted.", methodName);
            return true;
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        if (batchSize <= 0) {
            throw new IllegalArgumentException(methodName + ": batchSize must be greater than 0 but was " + batchSize + ".");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        Descriptor.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

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
                            boolean success = Descriptor.setDescriptorsForSingleMoleculeSmilesString(
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
                            Descriptor.LOGGER.log(
                                    Level.WARNING,
                                    String.format("Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization: Exception in batch %d, molecule index: %d", batchIndex, i),
                                    exception
                            );
                        }
                    }
                });

            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, nanPositionsList)) {
                        hasNaN.set(true);
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            Descriptor.LOGGER.log(
                    Level.WARNING,
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
        if (!Descriptor.validateDescriptors(descriptors, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given descriptor array is empty, calculation aborted.", methodName);
            return true;
        }
        if (!Descriptor.validateAtomContainerArray(atomContainerArray, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given atom container array is empty, calculation aborted.", methodName);
            return true;
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        Descriptor.validateMatrix(matrix, descriptors, atomContainerArray, startIndex, methodName);

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
                                boolean success = Descriptor.setDescriptorsForSingleMolecule(
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
                                Descriptor.LOGGER.log(
                                        Level.WARNING,
                                        String.format("Descriptor.setDescriptorsForMoleculesByMoleculeParallelization: One descriptor calculation caused an exception, molecule index: %d.", i),
                                        exception
                                );
                            }
                        }
                );
            } else {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    if (!Descriptor.setDescriptorsForSingleMolecule(descriptors, atomContainerArray[i], matrix[i], startIndices, i, nanPositionsList)) {
                        hasNaN.set(true);
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            Descriptor.LOGGER.log(
                    Level.WARNING,
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
        if (!Descriptor.validateDescriptors(descriptors, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given descriptor array is empty, calculation aborted.", methodName);
            return true;
        }
        if (!Descriptor.validateSmilesStringArray(moleculeSmilesStringArray, methodName)) {
            Descriptor.LOGGER.log(Level.WARNING, "{0} : Given SMILES string array is empty, calculation aborted.", methodName);
            return true;
        }
        if (nanPositionsList == null) {
            throw new NullPointerException(methodName + ": nanPositionsList is null.");
        }
        // throws NullPointerException or IllegalArgumentException if the matrix or on eof its rows is null or its dimensions are invalid
        Descriptor.validateMatrix(matrix, descriptors, moleculeSmilesStringArray, startIndex, methodName);

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
                                boolean success = Descriptor.setDescriptorsForSingleMoleculeSmilesString(
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
                                Descriptor.LOGGER.log(
                                        Level.WARNING,
                                        String.format("Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization: Exception in molecule index: %d", i),
                                        exception
                                );
                            }

                        }
                );
            } else {
                for (int i = 0; i < numberOfMolecules; i++) {
                    if (!Descriptor.setDescriptorsForSingleMoleculeSmilesString(descriptors, moleculeSmilesStringArray[i], matrix[i], startIndices, i, electronDonationModel, nanPositionsList)) {
                        hasNaN.set(true);
                    }
                }
            }
        } catch (InterruptedException exception) {
            throw exception;
        } catch (Exception exception) {
            Descriptor.LOGGER.log(
                    Level.WARNING,
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
                moleculeWithExplicitHydrogens = Descriptor.createMoleculeWithExplicitHydrogens(atomContainer);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : atomContainer;
            if (!Descriptor.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
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
        IAtomContainer molecule = Descriptor.SMILES_PARSER.parseSmiles(moleculeSmilesString);
        if (electronDonationModel == null) {
            Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        } else {
            Descriptor.setAromaticity(molecule, electronDonationModel);
        }

        IAtomContainer moleculeWithExplicitHydrogens = null;
        for (Descriptor descriptor : descriptors) {
            if (descriptor.needsExplicitHydrogens()){
                moleculeWithExplicitHydrogens = Descriptor.createMoleculeWithExplicitHydrogens(molecule);
                break;
            }
        }
        boolean isSuccessful = true;
        for (int i = 0; i < descriptors.length; i++) {
            IAtomContainer moleculeToUse = descriptors[i].needsExplicitHydrogens() ? moleculeWithExplicitHydrogens : molecule;
            if (!Descriptor.setDescriptor(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
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
            boolean success = descriptor.calculate(atomContainer, vector, startIndex);
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
            return !Descriptor.checkAndTrackNaNValues(vector, startIndex, numComponents, moleculeIndex, nanPositionsList);
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
            Descriptor.LOGGER.log(
                    Level.WARNING,
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
            return Descriptor.copyMolecule(molecule);
        }
        try {
            // Create a deep copy of the molecule first
            IAtomContainer moleculeCopy = Descriptor.copyMolecule(molecule);
            // Add explicit hydrogen atoms
            AtomContainerManipulator.normalizeHydrogens(moleculeCopy, HydrogenState.Explicit);
            return moleculeCopy;
        } catch (Exception exception) {
            throw new CloneNotSupportedException("Could not create molecule with explicit hydrogens: " + exception.getMessage());
        }
    }
}
