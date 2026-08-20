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


import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
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
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.geometry.GeometryUtil;
import org.openscience.cdk.io.HINReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WHIMDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WeightedPathDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.WienerNumbersDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.ZagrebIndexDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.qsar.result.DoubleResult;
import org.openscience.cdk.qsar.result.IntegerArrayResult;
import org.openscience.cdk.qsar.result.IntegerResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import java.io.InputStream;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.atomic.AtomicBoolean;
import java.util.stream.IntStream;

//TODO test descriptors for how they handle empty molecules and empty SMILES strings.

/**
 * Test class for Descriptor class. This class first tests all the descriptors included in the {@link Descriptor} class
 * individually, i.e. whether they produce the expected results for some example molecules (in most cases taken from the
 * respective CDK descriptor test class) in different parallelization settings.
 * Note: For adding tests of a new descriptor goto "Add new descriptor tests here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorTest {
    /**
     * Test method for descriptor MOLECULAR_WEIGHT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MOLECULAR_WEIGHT() throws Exception {
        String smiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 60.05f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor WIENER_NUMBER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_WIENER_NUMBER() throws Exception {
        String smiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.WIENER_NUMBER};
        boolean isParallelCalculation = false;
        float expectedWienerPath = 9.0f; // 1(C1C2)+2(C1O1)+2(C1O2)+1(C2O1)+1(C2O2)+2(O1O2) = 9
        float expectedWienerPolarity = 0.0f; // there are no atoms that are 3 bonds apart

        Assertions.assertEquals(2, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedWienerPath, matrix[0][0]);
        Assertions.assertEquals(expectedWienerPolarity, matrix[0][1]);

        matrix = new float[][]
                {
                        {0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedWienerPath, matrix[0][0]);
        Assertions.assertEquals(expectedWienerPolarity, matrix[0][1]);
    }

    /**
     * Test method for descriptor ATOM_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ATOM_COUNT() throws Exception {
        String smiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ATOM_COUNT};
        boolean isParallelCalculation = false;
        float expected = 8.0f; // Acetic acid has 8 atoms (implicit Hs included)

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };

        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor ATOM_COUNT_HEAVY.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ATOM_COUNT_HEAVY() throws Exception {
        String smiles = "CC(=O)O"; //Acetic Acid CID: 176
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ATOM_COUNT_HEAVY};
        boolean isParallelCalculation = false;
        float expected = 4.0f; // Acetic acid has 4 heavy atoms

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };

        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for the organic subset of individual atom counts (C, H, N, O, S, P, F, Br, Cl, I) in a complex molecule.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ATOM_COUNT_ORGANIC_SUBSET() throws Exception {
        // Cobalamin CID: 74413906
        String smiles = "CC1=CC2=C(C=C1C)N(C=N2)C3C(C(C(O3)CO)OP(=O)([O-])OC(C)CNC(=O)CCC4(C(C5C6(C(C(C(=N6)C(=C7C(C(C(=N7)C=C8C(C(C(=N8)C(=C4[N-]5)C)CCC(=O)N)(C)C)CCC(=O)N)(C)CC(=O)N)C)CCC(=O)N)(C)CC(=O)N)C)CC(=O)N)C)O.[Co+3]";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.ATOM_COUNT_C,
                Descriptor.ATOM_COUNT_H,
                Descriptor.ATOM_COUNT_N,
                Descriptor.ATOM_COUNT_O,
                Descriptor.ATOM_COUNT_S,
                Descriptor.ATOM_COUNT_P,
                Descriptor.ATOM_COUNT_F,
                Descriptor.ATOM_COUNT_BR,
                Descriptor.ATOM_COUNT_CL,
                Descriptor.ATOM_COUNT_I
        };
        boolean isParallelCalculation = false;

        Assertions.assertEquals(10, Descriptor.getNumberOfComponents(descriptors));

        float[] expected = new float[]{
                62, // C count
                88, // H count
                13, // N count
                14, // O count
                0, // S count
                1, // P count
                0, // F count
                0, // Br count
                0, // Cl count
                0 // I count
        };

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f,}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0]);
    }

    /**
     * Test method for descriptor H_BOND_ACCEPTOR_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_H_BOND_ACCEPTOR_COUNT() throws Exception {
        String smiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        String smiles2 = "O=N(=O)c1cccc2cn[nH]c12"; // 7-Nitroindole CID: 1893
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.H_BOND_ACCEPTOR_COUNT};
        boolean isParallelCalculation = false;
        float expectedMol1 = 2f;
        float expectedMol2 = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                // Mdl and PiBonds lead to a calculation of 2 rather than 1 H-bond acceptor which is correct
                // because of the different handling of aromaticity in the second molecule but this behavior
                // is not tested here
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                //Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles
                //Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);

            matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);
        }
    }

    /**
     * Test method for descriptor H_BOND_DONOR_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_H_BOND_DONOR_COUNT() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        String smiles1 = "CC(=O)O"; //Acetic Acid CID: 176
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        String smiles2 = "Oc1ccccc1"; // Phenol CID: 996
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.H_BOND_DONOR_COUNT};
        boolean isParallelCalculation = false;
        float expectedMol1 = 1f; // Acetic acid has 1 hydrogen bond donor
        float expectedMol2 = 1f; // Phenol has 1 hydrogen bond donor

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);

            matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);
        }
    }

    /**
     * Test method for descriptor TPSA.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_TPSA() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // 3-(dimethylamino)-3-(methyleneamino)propanenitrile (not in PubChem)
        String smiles1 = "C=NC(CC#N)N(C)C";
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        // 1-Nitropropane CID: 7903
        String smiles2 = "CCCN(=O)=O";
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.TPSA};
        boolean isParallelCalculation = false;

        double epsilon = 0.01; //tolerance range
        float expectedMol1 = 39.39f;
        float expectedMol2 = 45.82f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0], epsilon);
            Assertions.assertEquals(expectedMol2, matrix[1][0], epsilon);

            matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0], epsilon);
            Assertions.assertEquals(expectedMol2, matrix[1][0], epsilon);
        }
    }

    /**
     * Test method for descriptor LARGEST_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LARGEST_CHAIN() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // 1-Phenylbutadiene CID: 137048
        String smiles1 = "C=CC=Cc1ccccc1";
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        // 4-(4-(penta-2,4-dien-1-yl)benzyl)-3-vinylpyridine (not in PubChem)
        String smiles2 = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2";
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.LARGEST_CHAIN};
        boolean isParallelCalculation = false;
        float expectedMol1 = 4f;
        float expectedMol2 = 5f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        //TODO: does this descriptor really need aromaticity info?
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);

            matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);
        }
    }

    /**
     * Test method for descriptor LONGEST_ALIPHATIC_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LONGEST_ALIPHATIC_CHAIN() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Butylbenzene CID: 7705
        String smiles1 = "CCCCc1ccccc1";
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        // 4-(4-tert-butylphenoxy)-N-(1,3-thiazol-2-yl)butanamide CID: 1565007
        String smiles2 = "CC(C)(C)c2ccc(OCCCC(=O)Nc1nccs1)cc2";
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.LONGEST_ALIPHATIC_CHAIN};
        boolean isParallelCalculation = false;
        float expectedMol1 = 4f;
        float expectedMol2 = 4f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model to the first molecule
            Descriptor.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            Descriptor.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);

            matrix = new float[][]
                    {
                            {0f}, {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expectedMol1, matrix[0][0]);
            Assertions.assertEquals(expectedMol2, matrix[1][0]);
        }
    }

    /**
     * Test method for descriptor MANNHOLD_LOGP.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MANNHOLD_LOGP() throws Exception {
        String smiles = "C"; // Methane CID: 297
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MANNHOLD_LOGP};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; //tolerance range
        float expected = 1.57f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor BCUT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BCUT() throws Exception {
        String smiles = "CC(=O)N"; // Acetamide CID: 178
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.BCUT};
        boolean isParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range

        float[] expected = new float[] {11.881587f, 16.005958f, -0.381844f, 0.325509f, 3.374638f, 5.033583f};

        Assertions.assertEquals(6, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor BOND_COUNT_ALL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BOND_COUNT_ALL() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // Ethanol CID: 702
        String smiles1 = "CCO";
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        // Allene CID: 10037
        String smiles2 = "C=C=C";
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);
        // Acetonitrile CID: 6342
        String smiles3 = "CC#N";
        IAtomContainer molecule3 = smilesParser.parseSmiles(smiles3);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2, molecule3};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.BOND_COUNT_ALL};
        boolean isParallelCalculation = false;
        //same total bond count for all three molecules
        float expected = 2f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}, {0f}, {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
        Assertions.assertEquals(expected, matrix[1][0]);
        Assertions.assertEquals(expected, matrix[2][0]);

        matrix = new float[][]
                {
                        {0f}, {0f}, {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
        Assertions.assertEquals(expected, matrix[1][0]);
        Assertions.assertEquals(expected, matrix[2][0]);
    }

    /**
     * Test method for descriptors BOND_COUNT_SINGLE, BOND_COUNT_DOUBLE, and BOND_COUNT_TRIPLE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BOND_COUNT_SPECIFIED() throws Exception {
        String smiles = "C=CC#N"; // Acrylonitrile CID: 7855
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.BOND_COUNT_SINGLE,
                Descriptor.BOND_COUNT_DOUBLE,
                Descriptor.BOND_COUNT_TRIPLE
        };
        boolean isParallelCalculation = false;

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(descriptors));

        float expectedSingleBonds = 1f;
        float expectedDoubleBonds = 1f;
        float expectedTripleBonds = 1f;

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedSingleBonds, matrix[0][0]);
        Assertions.assertEquals(expectedDoubleBonds, matrix[0][1]);
        Assertions.assertEquals(expectedTripleBonds, matrix[0][2]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedSingleBonds, matrix[0][0]);
        Assertions.assertEquals(expectedDoubleBonds, matrix[0][1]);
        Assertions.assertEquals(expectedTripleBonds, matrix[0][2]);
    }

    /**
     * Test method for descriptor B_POL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_B_POL() throws Exception {
        String smiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.B_POL};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 7.517242f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor RULE_OF_FIVE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_RULE_OF_FIVE() throws Exception {
        // 3-(4-(2-(1-ethoxybutoxy)-1-(naphthalen-1-yl)ethyl)-2-((hexahydropyrimidin-5-yl)methyl)cyclohexyl)propan-1-ol (not in PubChem)
        String smiles = "CCCC(OCC)OCC(c1cccc2ccccc12)C4CCC(CCCO)C(CC3CNCNC3)C4";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.RULE_OF_FIVE};
        boolean isParallelCalculation = false;
        float expected = 3f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);

            matrix = new float[][]
                    {
                            {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);
        }
    }

    /**
     * Test method for descriptor AROMATIC_ATOMS_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AROMATIC_ATOMS_COUNT() throws Exception {
        String smiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AROMATIC_ATOMS_COUNT};
        boolean isParallelCalculation = false;
        float expected = 6f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        // Array of all available electron donation models
        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);

            matrix = new float[][]
                    {
                            {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);
        }
    }

    /**
     * Test method for descriptor AROMATIC_BONDS_COUNT with benzene.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AROMATIC_BONDS_COUNT() throws Exception {
        String smiles = "c1ccccc1"; // Benzene CID: 241
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AROMATIC_BONDS_COUNT};
        boolean isParallelCalculation = false;
        float expected = 6f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);

            matrix = new float[][]
                    {
                            {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0]);
        }
    }

    /**
     * Test method for descriptor ROTATABLE_BONDS_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ROTATABLE_BONDS_COUNT() throws Exception {
        // N-ethyl-3-methylbutanamide CID: 528605
        String smiles = "CCNC(=O)CC(C)C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ROTATABLE_BONDS_COUNT};
        boolean isParallelCalculation = false;
        float expected = 4f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor FMF.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FMF() throws Exception {
        // Clenbuterol CID: 2783
        String smiles = "Clc1cc(cc(Cl)c1N)C(O)CNC(C)(C)C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.FMF};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 0.353f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor FRACTIONAL_CSP3.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRACTIONAL_CSP3() throws Exception {
        // 2,6-Dimethylpyridine CID: 7937
        String smiles = "CC1=CC=CC(C)=N1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.FRACTIONAL_CSP3};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 0.29f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor HYBRIDIZATION_RATIO.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_HYBRIDIZATION_RATIO() throws Exception {
        String smiles = "CCC"; // Propane CID: 6334
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.HYBRIDIZATION_RATIO};
        boolean isParallelCalculation = false;
        float expected = 1.00f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor KAPPA_SHAPE_INDICES.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_KAPPA_SHAPE_INDICES() throws Exception {
        String smiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.KAPPA_SHAPE_INDICES};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] expected = new float[] {5.0f, 2.25f, 4.0f};

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor PETITJEAN_NUMBER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_PETITJEAN_NUMBER() throws Exception {
        String smiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.PETITJEAN_NUMBER};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 0.33333334f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor SPIRO_ATOM_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_SPIRO_ATOM_COUNT() throws Exception {
        // 3'H-spiro[cyclohexane-1,2'-naphthalene] (not in PubChem)
        String smiles = "C1CCC2(CC1)CC=C1C=CC=CC1=C2";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.SPIRO_ATOM_COUNT};
        boolean isParallelCalculation = false;
        float expected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor V_ADJ_MAT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_V_ADJ_MAT() throws Exception {
        String smiles = "C1CCC2CCCCC2C1"; // Decalin CID: 7044
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.V_ADJ_MAT};
        boolean isParallelCalculation = false;
        double epsilon = 0.001; // tolerance range
        float expected = 4.459f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor WEIGHTED_PATH.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_WEIGHTED_PATH() throws Exception {
        String smiles = "CCCC"; // Butane CID: 7843
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.WEIGHTED_PATH};
        boolean isParallelCalculation = false;
        double epsilon = 0.00001; // tolerance range

        float[] expected = new float[] {6.87132f, 1.71783f, 0.0f, 0.0f, 0.0f};

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor ZAGREB_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ZAGREB_INDEX() throws Exception {
        String smiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ZAGREB_INDEX};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range
        float expected = 16f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor CARBON_TYPES.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CARBON_TYPES() throws Exception {
        String smiles = "CCCC"; // Butane CID: 7843
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CARBON_TYPES};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(9, Descriptor.getNumberOfComponents(descriptors));

        float[] expected = new float[] {0f, 0f, 0f, 0f, 0f, 2f, 2f, 0f, 0f};

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0]);
    }

    /**
     * Test method for descriptor A_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_A_LOG_P() throws Exception {
        String smiles = "CCCCl"; // 1-Chloropropane CID: 10899
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[] {molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[] {Descriptor.A_LOG_P};
        boolean isParallelCalculation = false;
        double epsilon = 0.001; //tolerance range

        float[] expected = new float[] {1.719f, 2.955f, 20.584f};

        Assertions.assertEquals(3, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor X_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_X_LOG_P() throws Exception {
        String smiles = "O=C(O)C(N)CCCN"; // 2,5-Diaminopentanoic Acid CID: 389
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.X_LOG_P};
        boolean isParallelCalculation = false;
        double epsilon = 0.1; // tolerance range
        float expected = -3.30f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[][]
                    {
                            {0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0], epsilon);

            matrix = new float[][]
                    {
                            {0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected, matrix[0][0], epsilon);
        }
    }

    /**
     * Test method for descriptor JP_LOG_P.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_JP_LOG_P() throws Exception {
        String smiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.JP_LOG_P};
        boolean isParallelCalculation = false;
        double epsilon = 0.1; // tolerance range
        float expected = 0.3f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor A_POL.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_A_POL() throws Exception {
        String smiles = "O=C(O)CC"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.A_POL};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 10.88f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f} //TODO: BTW, does the matrix necessarily have to be initialised with zeros to use the library?
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_CHARGE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_CHARGE() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_CHARGE};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO: I know, this is taken like that from the CDK tests; but could we maybe test for the actual values instead of != 0?
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_MASS.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_MASS() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_MASS};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO see above
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_POLARIZABILITY.
     * No validated result because the descriptor itself is not validated in the CDK.
     * Result can be printed to see if descriptor calculates values.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_POLARIZABILITY() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_POLARIZABILITY};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO: we need to find a better solution for this, there should be no print-outs; maybe there are example values in the original publication?
        System.out.println("New calculation results:");
        System.out.println(Arrays.toString(matrix[0]));

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        System.out.println("Descriptor parallelization results:");
        System.out.println(Arrays.toString(matrix[0]));
    }

    /**
     * Test method for descriptor FRAGMENT_COMPLEXITY.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRAGMENT_COMPLEXITY() throws Exception {
        String smiles = "c1ccc(CCc2ccccc2)cc1"; // 1,2-Diphenylethane CID: 7647
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.FRAGMENT_COMPLEXITY};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 659.00f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor CHI_CHAIN.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_CHAIN() throws Exception {
        String smiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CHI_CHAIN};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] expected = new float[] {0.2887f, 0.2887f, 0.0000f, 0.0000f, 0.0000f, 0.1667f, 0.1667f, 0.0000f, 0.0000f};

        Assertions.assertEquals(10, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_CLUSTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_CLUSTER() throws Exception {
        String smiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CHI_CLUSTER};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] expected = new float[]{0.2887f, 0.0000f, 0.0000f, 0.0000f, 0.1667f, 0.0000f, 0.0000f, 0.0000f};

        Assertions.assertEquals(8, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_PATH_CLUSTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_PATH_CLUSTER() throws Exception {
        String smiles = "C1=C(Cl)C=CC=C1(Cl)"; // 1,3-Dichlorobenzene CID: 10943
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CHI_PATH_CLUSTER};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] expected = new float[]{0.7416f, 1.0934f, 1.0202f, 0.4072f, 0.5585f, 0.4376f};

        Assertions.assertEquals(6, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor CHI_PATH.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CHI_PATH() throws Exception {
        String smiles = "CC1OC1"; // Propylene oxide CID: 6378
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CHI_PATH};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        float[] expected = new float[] {2.9916f, 1.8938f, 1.6825f, 0.5773f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, 2.6927f, 1.5099f, 1.1439f};

        Assertions.assertEquals(16, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        for (int i = 0; i < expected.length; i++) {
            //cannot use Assertions.assertArrayEquals() because there is no epsilon
            Assertions.assertEquals(expected[i], matrix[0][i], epsilon);
        }
    }

    /**
     * Test method for descriptor FRACTIONAL_PSA.
     * Expected results were calculated by TPSADescriptor / MolecularWeightDescriptor.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_FRACTIONAL_PSA() throws Exception {
        String smiles = "O=C(O)c1ccncc1"; // Isonicotinic Acid CID: 5922
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.FRACTIONAL_PSA};
        boolean isParallelCalculation = false;
        double epsilon = 0.001; // tolerance range
        float expected = 0.4077f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Test method for descriptor LARGEST_PI_SYSTEM.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LARGEST_PI_SYSTEM() throws Exception {
        //4-(4-(penta-2,4-dien-1-yl)benzyl)-3-vinylpyridine (not in PubChem)
        String smiles = "C=CC=CCc2ccc(Cc1ccncc1C=C)cc2";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.LARGEST_PI_SYSTEM};
        boolean isParallelCalculation = false;
        float expected = 8f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

    }

    /**
     * Test method for descriptor SMALL_RING.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_SMALL_RING() throws Exception {
        // 5,14-Pentacenedione CID: 10686237
        String smiles = "O=C1c2ccccc2C(=O)c2cc3cc4ccccc4cc3cc21";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);

        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.SMALL_RING};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(11, Descriptor.getNumberOfComponents(descriptors));

        //numbers indicate the tested descriptor value, not all values are tested here, only the first 4 positions
        float expected0 = 5f;
        float expected1 = 5f;
        float expected2 = 1f;
        float expected3 = 1f;

        ElectronDonation[] models = {
                Aromaticity.Model.Daylight,
                Aromaticity.Model.CDK_2x,
                Aromaticity.Model.CDK_1x,
                Aromaticity.Model.CDK_AtomTypes,
                Aromaticity.Model.Mdl,
                Aromaticity.Model.OpenSmiles,
                Aromaticity.Model.PiBonds
        };

        for (ElectronDonation model : models) {
            // Apply aromaticity with the current model
            Descriptor.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected0, matrix[0][0]);
            Assertions.assertEquals(expected1, matrix[0][1]);
            Assertions.assertEquals(expected2, matrix[0][2]);
            Assertions.assertEquals(expected3, matrix[0][3]);

            matrix = new float[][]
                    {
                            {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                    };
            nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                            descriptors,
                            moleculesArray,
                            matrix,
                            startIndex,
                            isParallelCalculation,
                            nanPositions
                    )
            );
            Assertions.assertEquals(expected0, matrix[0][0]);
            Assertions.assertEquals(expected1, matrix[0][1]);
            Assertions.assertEquals(expected2, matrix[0][2]);
            Assertions.assertEquals(expected3, matrix[0][3]);
        }
    }

    /**
     * Test method for descriptor BASIC_GROUP_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_BASIC_GROUP_COUNT() throws Exception {
        String smiles = "NC"; // Methylamine CID: 6329
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.BASIC_GROUP_COUNT};
        boolean isParallelCalculation = false;
        float expected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor ACIDIC_GROUP_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ACIDIC_GROUP_COUNT() throws Exception {
        String smiles = "CC(=O)O"; // Acetic acid CID: 176
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ACIDIC_GROUP_COUNT};
        boolean isParallelCalculation = false;
        float expected = 1f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor AMINO_ACID_COUNT.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AMINO_ACID_COUNT() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        // L-threonyl-L-threonine CID: 11321969
        String smiles1 = "N[C@@]([H])([C@]([H])(O)C)C(=O)N[C@@]([H])([C@]([H])(O)C)C(=O)O";
        IAtomContainer molecule1 = smilesParser.parseSmiles(smiles1);
        // Glycylglycine CID: 11163
        String smiles2 = "C(C(=O)NCC(=O)O)N";
        IAtomContainer molecule2 = smilesParser.parseSmiles(smiles2);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AMINO_ACID_COUNT};
        boolean isParallelCalculation = false;
        //expected threonine count does not apply to the 2nd mol
        float expectedGlycineAndThreonineCount = 2f;

        Assertions.assertEquals(20, Descriptor.getNumberOfComponents(descriptors));


        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[0][8]); //TODO: is this the inherent behaviour of this descriptor, that all amino acid moieties are also detected as Glycine?
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[0][16]);
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[1][8]);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f},
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[0][8]);
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[0][16]);
        Assertions.assertEquals(expectedGlycineAndThreonineCount, matrix[1][8]);

    }

    /**
     * Test method for descriptor KIER_HALL_SMARTS.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_KIER_HALL_SMARTS() throws Exception {
        //5-(3-(aminomethyl)-5-(2-(methylamino)ethyl)phenyl)-6-hydroxy-1-(3-hydroxypropoxy)hexan-3-one (not in PubChem)
        String smiles = "c1c(CN)cc(CCNC)cc1C(CO)CC(=O)CCOCCCO";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.KIER_HALL_SMARTS};
        boolean isParallelCalculation = false;
        // numbers indicate the vector positions of the tested descriptor values; not all return values are tested here
        float expected33 = 2f;
        float expected34 = 1f;
        float expected35 = 1f;
        float expected20 = 1f;
        float expected23 = 1f;

        Assertions.assertEquals(79, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        new float[79]
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected33, matrix[0][33]);
        Assertions.assertEquals(expected34, matrix[0][34]);
        Assertions.assertEquals(expected35, matrix[0][35]);
        Assertions.assertEquals(expected20, matrix[0][20]);
        Assertions.assertEquals(expected23, matrix[0][23]);

        matrix = new float[][]
                {
                        new float[79]
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected33, matrix[0][33]);
        Assertions.assertEquals(expected34, matrix[0][34]);
        Assertions.assertEquals(expected35, matrix[0][35]);
        Assertions.assertEquals(expected20, matrix[0][20]);
        Assertions.assertEquals(expected23, matrix[0][23]);
    }

    /**
     * Test method for descriptor ECCENTRIC_CONNECTIVITY_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_ECCENTRIC_CONNECTIVITY_INDEX() throws Exception {
        // (1R,4S,5S,8S,9R,12S,13R)-1,5,9-trimethyl-11,14,15,16-tetraoxatetracyclo[10.3.1.04,13.08,13]hexadecan-10-one CID: 98047509
        String smiles = "C[C@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.ECCENTRIC_CONNECTIVITY_INDEX};
        boolean isParallelCalculation = false;
        float expected = 254f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0]);
    }

    /**
     * Test method for descriptor MDE.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MDE() throws Exception {
        //2,2-bis(methylperoxy)propan-1-ol (not in PubChem)
        String smiles = "COOC(C)(CO)OOC";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MDE};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001; // tolerance range

        //numbers indicate the position of the tested descriptor result; not all result values are tested here
        float expected10 = 0.0000f;
        float expected11 = 1.1547f;
        float expected12 = 2.9416f;

        Assertions.assertEquals(19, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected10, matrix[0][10], epsilon);
        Assertions.assertEquals(expected11, matrix[0][11], epsilon);
        Assertions.assertEquals(expected12, matrix[0][12], epsilon);

        matrix = new float[][]
                {
                        {0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f, 0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected10, matrix[0][10], epsilon);
        Assertions.assertEquals(expected11, matrix[0][11], epsilon);
        Assertions.assertEquals(expected12, matrix[0][12], epsilon);
    }

    /**
     * Test method for descriptor VABC.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_VABC() throws Exception {
        // Prilosec CID: 4594
        String smiles = "COc2ccc1[nH]c(nc1c2)S(=O)Cc3ncc(C)c(OC)c3C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.VABC};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 292.23f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[][]
                {
                        {0f}
                };
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);

        matrix = new float[][]
                {
                        {0f}
                };
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(expected, matrix[0][0], epsilon);
    }

    /**
     * Tests PubChem fingerprint descriptor.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_PUBCHEM_FINGERPRINT() throws Exception {
        // 1-Benzyl-2,4-diphenyl-6-(2-phenylethenyl)pyridin-1-ium CID: 3828524
        String smiles = "C1=CC=C(C=C1)C[N+]2=C(C=C(C=C2C=CC3=CC=CC=C3)C4=CC=CC=C4)C5=CC=CC=C5";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER};
        boolean isParallelCalculation = false;
        // Verify that expected bits are set (convert reference fingerprint to expected values)
        BitSet ref = PubchemFingerprinter.decode(
                "AAADceB+AAAAAAAAAAAAAAAAAAAAAAAAAAA8YMGCAAAAAAAB1AAAHAAAAAAADAjBHgQwgJMMEACgAyRiRACCgCAhAiAI2CA4ZJgIIOLAkZGEIAhggADIyAcQgMAOgAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=="
        );

        // Test descriptor component count
        Assertions.assertEquals(881, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][881];
        Arrays.fill(matrix[0], 0f);
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        for (int i = 0; i < 881; i++) {
            float expected = ref.get(i) ? 1.0f : 0.0f;
            Assertions.assertEquals(expected, matrix[0][i]);
        }

        matrix = new float[1][881];
        Arrays.fill(matrix[0], 0f);
        nanPositions = Collections.synchronizedList(new LinkedList<>());

        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        for (int i = 0; i < 881; i++) {
            float expected = ref.get(i) ? 1.0f : 0.0f;
            Assertions.assertEquals(expected, matrix[0][i]);
        }

    }

    /**
     * Test method for descriptor MACCS_FINGERPRINTER.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MACCS_FINGERPRINTER() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String smiles = "c1ccccc1CCc1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MACCS_FINGERPRINTER};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(166, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][166];
        Arrays.fill(matrix[0], 0f);
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO: is there no reference fingerprint?
        Assertions.assertEquals(1, matrix[0][124]);
        Assertions.assertEquals(0, matrix[0][165]);

        matrix = new float[1][166];
        Arrays.fill(matrix[0], 0f);
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(1, matrix[0][124]);
        Assertions.assertEquals(0, matrix[0][165]);
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_ECFP descriptors of different diameters.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_ECFP() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String smiles = "c1ccccc1CCc1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6
        };
        boolean isParallelCalculation = false;

        Assertions.assertEquals(1024 * 4, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][1024 * 4];
        Arrays.fill(matrix[0], 0f);
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO see above
        Assertions.assertEquals(0, nanPositions.size());

        matrix = new float[1][1024 * 4];
        Arrays.fill(matrix[0], 0f);
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0, nanPositions.size());
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_FCFP descriptors of different diameters.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_FCFP() throws Exception {
        //1,2-Diphenylethane CID: 7647
        String smiles = "c1ccccc1CCc1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6
        };
        boolean isParallelCalculation = false;

        Assertions.assertEquals(4 * 1024, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][4 * 1024];
        Arrays.fill(matrix[0], 0f);
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        //TODO see above
        Assertions.assertEquals(0, nanPositions.size());

        matrix = new float[1][4 * 1024];
        Arrays.fill(matrix[0], 0f);
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0, nanPositions.size());
    }

    /**
     * Test method for descriptor CPSA.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CPSA() throws Exception {
        String filename = "benzene.hin";
        InputStream ins = CPSADescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new HINReader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer molecule = (IAtomContainer) cList.get(0);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CPSA};
        boolean isParallelCalculation = false;
        double epsilon = 0.0001;

        Assertions.assertEquals(29, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][29];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0.0f, matrix[0][28], epsilon); // RPSA
        Assertions.assertEquals(1.0f, matrix[0][27], epsilon); // RHSA
        Assertions.assertEquals(0.0f, matrix[0][26], epsilon); // TPSA
        Assertions.assertEquals(231.66182f, matrix[0][25], epsilon); // THSA

        matrix = new float[1][29];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0.0f, matrix[0][28], epsilon);
        Assertions.assertEquals(1.0f, matrix[0][27], epsilon);
        Assertions.assertEquals(0.0f, matrix[0][26], epsilon);
        Assertions.assertEquals(231.66182f, matrix[0][25], epsilon);

        float[] singleResult = Descriptor.calculateDescriptor(Descriptor.CPSA, molecule);
        Assertions.assertEquals(0.0f, singleResult[28], epsilon);
        Assertions.assertEquals(1.0f, singleResult[27], epsilon);
        Assertions.assertEquals(0.0f, singleResult[26], epsilon);
        Assertions.assertEquals(231.66182f, singleResult[25], epsilon);
    }

    /**
     * Test method for descriptor GRAVITATIONAL_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_GRAVITATIONAL_INDEX() throws Exception {
        String filename = "gravindex.hin";
        InputStream ins = GravitationalIndexDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new HINReader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer molecule = (IAtomContainer) cList.get(0);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.GRAVITATIONAL_INDEX};
        boolean isParallelCalculation = false;
        double epsilon = 0.00001;

        Assertions.assertEquals(9, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][9];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(1756.50607f, matrix[0][0], epsilon);
        Assertions.assertEquals(41.91069f, matrix[0][1], epsilon);
        Assertions.assertEquals(12.06563f, matrix[0][2], epsilon);
        Assertions.assertEquals(1976.64326f, matrix[0][3], epsilon);
        Assertions.assertEquals(44.45946f, matrix[0][4], epsilon);
        Assertions.assertEquals(12.54997f, matrix[0][5], epsilon);
        Assertions.assertEquals(4333.09737f, matrix[0][6], epsilon);
        Assertions.assertEquals(65.82627f, matrix[0][7], epsilon);
        Assertions.assertEquals(16.30295f, matrix[0][8], epsilon);

        matrix = new float[1][9];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(1756.50607f, matrix[0][0], epsilon);
        Assertions.assertEquals(41.91069f, matrix[0][1], epsilon);
        Assertions.assertEquals(12.06563f, matrix[0][2], epsilon);
        Assertions.assertEquals(1976.64326f, matrix[0][3], epsilon);
        Assertions.assertEquals(44.45946f, matrix[0][4], epsilon);
        Assertions.assertEquals(12.54997f, matrix[0][5], epsilon);
        Assertions.assertEquals(4333.09737f, matrix[0][6], epsilon);
        Assertions.assertEquals(65.82627f, matrix[0][7], epsilon);
        Assertions.assertEquals(16.30295f, matrix[0][8], epsilon);

        float[] singleResult = Descriptor.calculateDescriptor(Descriptor.GRAVITATIONAL_INDEX, molecule);
        Assertions.assertEquals(1756.50607f, singleResult[0], epsilon);
        Assertions.assertEquals(41.91069f, singleResult[1], epsilon);
        Assertions.assertEquals(12.06563f, singleResult[2], epsilon);
    }

    /**
     * Test method for descriptor MOMENT_OF_INERTIA.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MOMENT_OF_INERTIA() throws Exception {
        String filename = "gravindex.hin";
        InputStream ins = MomentOfInertiaDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new HINReader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer molecule = (IAtomContainer) cList.get(0);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOMENT_OF_INERTIA};
        boolean isParallelCalculation = false;
        double epsilon = 0.00001;

        Assertions.assertEquals(7, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][7];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(1820.692519f, matrix[0][0], epsilon);
        Assertions.assertEquals(1274.532522f, matrix[0][1], epsilon);
        Assertions.assertEquals(979.210423f, matrix[0][2], epsilon);
        Assertions.assertEquals(1.428517f, matrix[0][3], epsilon);
        Assertions.assertEquals(1.859347f, matrix[0][4], epsilon);
        Assertions.assertEquals(1.301592f, matrix[0][5], epsilon);
        Assertions.assertEquals(5.411195f, matrix[0][6], epsilon);

        matrix = new float[1][7];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(1820.692519f, matrix[0][0], epsilon);
        Assertions.assertEquals(1274.532522f, matrix[0][1], epsilon);
        Assertions.assertEquals(979.210423f, matrix[0][2], epsilon);
        Assertions.assertEquals(1.428517f, matrix[0][3], epsilon);
        Assertions.assertEquals(1.859347f, matrix[0][4], epsilon);
        Assertions.assertEquals(1.301592f, matrix[0][5], epsilon);
        Assertions.assertEquals(5.411195f, matrix[0][6], epsilon);

        float[] singleResult = Descriptor.calculateDescriptor(Descriptor.MOMENT_OF_INERTIA, molecule);
        Assertions.assertEquals(1820.692519f, singleResult[0], epsilon);
        Assertions.assertEquals(1274.532522f, singleResult[1], epsilon);
        Assertions.assertEquals(979.210423f, singleResult[2], epsilon);
    }

    /**
     * Test method for descriptor LENGTH_OVER_BREADTH.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_LENGTH_OVER_BREADTH() throws Exception {
        String filename = "lobtest.sdf";
        InputStream ins = LengthOverBreadthDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer molecule = (IAtomContainer) cList.get(0); // Cholesterol
        Isotopes.getInstance().configureAtoms(molecule);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.LENGTH_OVER_BREADTH};
        boolean isParallelCalculation = false;
        double epsilon = 0.001;

        Assertions.assertEquals(2, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][2];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(3.5029f, matrix[0][0], epsilon); // LOBMAX
        Assertions.assertEquals(3.5029f, matrix[0][1], epsilon); // LOBMIN

        matrix = new float[1][2];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(3.5029f, matrix[0][0], epsilon);
        Assertions.assertEquals(3.5029f, matrix[0][1], epsilon);

        float[] singleResult = Descriptor.calculateDescriptor(Descriptor.LENGTH_OVER_BREADTH, molecule);
        Assertions.assertEquals(3.5029f, singleResult[0], epsilon);
        Assertions.assertEquals(3.5029f, singleResult[1], epsilon);
    }

    /**
     * Test method for descriptor PETITJEAN_SHAPE_INDEX.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_PETITJEAN_SHAPE_INDEX() throws Exception {
        String filename = "petitejean.sdf";
        InputStream ins = PetitjeanShapeIndexDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new MDLV2000Reader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer molecule = (IAtomContainer) cList.get(0); // nbutane
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.PETITJEAN_SHAPE_INDEX};
        boolean isParallelCalculation = false;
        double epsilon = 0.000001;

        Assertions.assertEquals(2, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[1][2];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0.5f, matrix[0][0], epsilon); // topoShape
        Assertions.assertEquals(0.606477f, matrix[0][1], epsilon); // geomShape

        matrix = new float[1][2];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertEquals(0.5f, matrix[0][0], epsilon);
        Assertions.assertEquals(0.606477f, matrix[0][1], epsilon);

        float[] singleResult = Descriptor.calculateDescriptor(Descriptor.PETITJEAN_SHAPE_INDEX, molecule);
        Assertions.assertEquals(0.5f, singleResult[0], epsilon);
        Assertions.assertEquals(0.606477f, singleResult[1], epsilon);
    }

    // Add new descriptor tests here!

    /**
     * Tests parallelization. TODO: this test can split up into 5. And please add a bit more doc.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    //TODO for CDK integration: tag @Tag("SlowTest") needs to be added here
    void test_Parallelization() throws Exception {
        String smiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int numberOfMolecules = 1000;
        IAtomContainer[] moleculesArray = new IAtomContainer[numberOfMolecules];
        String[] moleculeStringsArray = new String[numberOfMolecules];

        //fill the arrays with 1000 (independent!) instances of propionic acid
        for (int i = 0; i < numberOfMolecules; i++) {
            IAtomContainer molecule = smilesParser.parseSmiles(smiles);
            Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
            moleculesArray[i] = molecule;
            moleculeStringsArray[i] = smiles;
        }
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);

        //first, calculate the results sequentially
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        //TODO: I again have the question here whether the matrix needs to be initialized with zeros or not; if not, we have to remove all my fill() calls above again and you should not initialize the arrays in the test methods above
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        //second, test the parallel calculation that uses a new descriptor instance for every calculation
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where single molecules are distributed onto threads
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where batches of molecules are distributed onto threads
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        int batchSize = 100;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where batches of SMILES strings are distributed onto threads
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixSequential,
                        startIndex,
                        batchSize,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixParallel,
                        startIndex,
                        batchSize,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        //override the results to now test the parallelization where single SMILES strings are distributed onto threads
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixSequential,
                        startIndex,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixParallel,
                        startIndex,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        //compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests integrity. TODO: please add a bit more doc.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    //TODO for CDK integration: tag @Tag("SlowTest") needs to be added here
    void test_Integrity() throws Exception {
        String smiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        int numberOfMolecules = 1000;
        IAtomContainer[] moleculesArray = new IAtomContainer[numberOfMolecules];
        String[] moleculeStringsArray = new String[numberOfMolecules];

        //fill the arrays with 1000 (independent!) instances of propionic acid
        for (int i = 0; i < numberOfMolecules; i++) {
            IAtomContainer molecule = smilesParser.parseSmiles(smiles);
            Descriptor.setAromaticity(molecule, Aromaticity.Model.Daylight);
            moleculesArray[i] = molecule;
            moleculeStringsArray[i] = smiles;
        }
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        boolean isParallelCalculation = false;
        int batchSize = 100;

        float[][] matrix1 = new float[numberOfMolecules][numberOfComponents];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix1,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrix2 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix2,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrix3 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrix3,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrix4 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrix4,
                        startIndex,
                        batchSize,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );
        float[][] matrix5 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesBySmilesStringParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrix5,
                        startIndex,
                        null, //use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );

        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrix1[i][j], matrix2[i][j]);
                Assertions.assertEquals(matrix1[i][j], matrix3[i][j]);
                Assertions.assertEquals(matrix1[i][j], matrix4[i][j]);
                Assertions.assertEquals(matrix1[i][j], matrix5[i][j]);
            }
        }
    }

    /**
     * Tests parallelization for 3D coordinates descriptors across calculation methods.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_Parallelization_3D() throws Exception {
        String filename = "gravindex.hin";
        InputStream ins = GravitationalIndexDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new HINReader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer baseMolecule = (IAtomContainer) cList.get(0);
        Isotopes.getInstance().configureAtoms(baseMolecule);

        int numberOfMolecules = 500;
        IAtomContainer[] moleculesArray = new IAtomContainer[numberOfMolecules];

        for (int i = 0; i < numberOfMolecules; i++) {
            moleculesArray[i] = baseMolecule.clone();
        }
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CPSA,
                Descriptor.GRAVITATIONAL_INDEX,
                Descriptor.MOMENT_OF_INERTIA,
                Descriptor.LENGTH_OVER_BREADTH,
                Descriptor.PETITJEAN_SHAPE_INDEX
        };
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);

        // First, calculate the results sequentially using setDescriptorsForMoleculesByMoleculeParallelizationNew
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        // Second, test the parallel calculation with new descriptor instances
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        // Compare results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        // Test parallelization by single molecule parallelization
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }

        // Test parallelization by batch parallelization
        matrixSequential = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = false;
        int batchSize = 10;
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositions
                )
        );

        matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests integrity for 3D coordinates descriptors across calculation methods.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_Integrity_3D() throws Exception {
        String filename = "gravindex.hin";
        InputStream ins = GravitationalIndexDescriptor.class.getResourceAsStream(filename);
        ISimpleChemObjectReader reader = new HINReader(ins);
        ChemFile content = (ChemFile) reader.read((ChemObject) new ChemFile());
        List<?> cList = ChemFileManipulator.getAllAtomContainers(content);
        IAtomContainer baseMolecule = (IAtomContainer) cList.get(0);
        Isotopes.getInstance().configureAtoms(baseMolecule);

        int numberOfMolecules = 500;
        IAtomContainer[] moleculesArray = new IAtomContainer[numberOfMolecules];

        for (int i = 0; i < numberOfMolecules; i++) {
            moleculesArray[i] = baseMolecule.clone();
        }
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CPSA,
                Descriptor.GRAVITATIONAL_INDEX,
                Descriptor.MOMENT_OF_INERTIA,
                Descriptor.LENGTH_OVER_BREADTH,
                Descriptor.PETITJEAN_SHAPE_INDEX
        };
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        boolean isParallelCalculation = false;
        int batchSize = 10;

        float[][] matrix1 = new float[numberOfMolecules][numberOfComponents];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew(
                        descriptors,
                        moleculesArray,
                        matrix1,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrix2 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix2,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrix3 = new float[numberOfMolecules][numberOfComponents];
        nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                Descriptor.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrix3,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositions
                )
        );

        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrix1[i][j], matrix2[i][j]);
                Assertions.assertEquals(matrix1[i][j], matrix3[i][j]);
            }
        }
    }

    /**
     * Test method for createMoleculeWithExplicitHydrogens.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testCreateMoleculeWithExplicitHydrogens() throws Exception {
        // Create a simple molecule (methane) with implicit hydrogen atoms
        String methaneSmiles = "C";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer methaneImplicit = smilesParser.parseSmiles(methaneSmiles);

        // Ensure the input is correct (1 atom, no bonds to H)
        Assertions.assertEquals(1, methaneImplicit.getAtomCount(),
                "Original methane should only have 1 heavy/explicit atom at this point");
        Assertions.assertEquals(0, methaneImplicit.getBondCount(),
                "Original methane should have no bonds");

        // Convert to molecule with explicit hydrogen atoms
        IAtomContainer methaneExplicit = Descriptor.createMoleculeWithExplicitHydrogens(methaneImplicit);

        // Check that the new molecule has the expected number of atoms (C + 4H = 5)
        Assertions.assertEquals(5, methaneExplicit.getAtomCount(),
                "Methane with explicit H should have 5 atoms (C + 4H)");

        // Check that the new molecule has the expected number of bonds (4 C-H bonds)
        Assertions.assertEquals(4, methaneExplicit.getBondCount(),
                "Methane with explicit H should have 4 bonds");

        // Check count of hydrogens
        int hydrogenCount = 0;
        for (IAtom atom : methaneExplicit.atoms()) {
            if ("H".equals(atom.getSymbol())) {
                hydrogenCount++;
            }
        }
        Assertions.assertEquals(4, hydrogenCount,
                "There should be 4 explicit hydrogen atoms");

        // Test with a more complex molecule
        String ethanolSmiles = "CCO";
        IAtomContainer ethanolImplicit = smilesParser.parseSmiles(ethanolSmiles);

        // Check if input is correct
        Assertions.assertEquals(3, ethanolImplicit.getAtomCount(),
                "Ethanol should have 3 heavy/explicit atoms here (all Hs are implicit)");

        // Convert to molecule with explicit hydrogens
        IAtomContainer ethanolExplicit = Descriptor.createMoleculeWithExplicitHydrogens(ethanolImplicit);

        // Check atom count (C + C + O + 6H = 9)
        Assertions.assertEquals(9, ethanolExplicit.getAtomCount(),
                "Ethanol with explicit H should have 9 atoms");
    }

    /**
     * Test method for copyMolecule method.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void testCopyMolecule() throws Exception {
        // Create a test molecule (benzene)
        String benzeneSmiles = "c1ccccc1";
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer benzene = smilesParser.parseSmiles(benzeneSmiles);

        // Verify initial state
        Assertions.assertEquals(6, benzene.getAtomCount(), "Original benzene should have 6 atoms");
        Assertions.assertEquals(6, benzene.getBondCount(), "Original benzene should have 6 bonds");

        // Create a copy of the molecule
        IAtomContainer benzeneClone = Descriptor.copyMolecule(benzene);

        // Verify the copy has the same structure
        Assertions.assertEquals(benzene.getAtomCount(), benzeneClone.getAtomCount(),
                "Clone should have the same number of atoms as the original");
        Assertions.assertEquals(benzene.getBondCount(), benzeneClone.getBondCount(),
                "Clone should have the same number of bonds as the original");

        // Ensure they are separate objects (deep copy)
        Assertions.assertNotSame(benzene, benzeneClone, "Clone should be a different object instance");

        // Modify the original to confirm the clone is independent
        IAtom newAtom = benzene.getBuilder().newInstance(IAtom.class, "O");
        benzene.addAtom(newAtom);

        // Verify the clone remains unchanged
        Assertions.assertEquals(7, benzene.getAtomCount(), "Modified original should have 7 atoms");
        Assertions.assertEquals(6, benzeneClone.getAtomCount(), "Clone should still have 6 atoms");
    }

    @Test
    void copyPreservesGenericAtomAndBondProperties() throws Exception {
        IAtomContainer original = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        IAtom a = original.newAtom(6, 3);
        IAtom b = original.newAtom(6, 3);
        a.setProperty("unique.atom.index", 0);
        b.setProperty("unique.atom.index", 1);
        IBond bond = original.newBond(a, b, IBond.Order.SINGLE);
        bond.setProperty("custom.bond.flag", "keep");
        original.setProperty("mol.title", "ethane");

        IAtomContainer copy = Descriptor.copyMolecule(original);

        Assertions.assertEquals(0, (int) copy.getAtom(0).getProperty("unique.atom.index"));
        Assertions.assertEquals(1, (int) copy.getAtom(1).getProperty("unique.atom.index"));
        Assertions.assertEquals("keep", copy.getBond(0).getProperty("custom.bond.flag"));
        Assertions.assertEquals("ethane", copy.getProperty("mol.title"));
        // Deep independence of the property map
        copy.getAtom(0).setProperty("unique.atom.index", 999);
        Assertions.assertEquals(0, (int) original.getAtom(0).getProperty("unique.atom.index"));
    }

    @Test
    void copyPreservesCoordinatesAsNewInstances() throws Exception {
        IAtomContainer original = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        IAtom a = original.newAtom(6, 4);
        a.setPoint2d(new Point2d(1.5, -2.0));
        a.setPoint3d(new Point3d(1.0, 2.0, 3.0));

        IAtomContainer copy = Descriptor.copyMolecule(original);
        IAtom c = copy.getAtom(0);

        Assertions.assertEquals(a.getPoint2d(), c.getPoint2d(), "2D coordinates must match");
        Assertions.assertEquals(a.getPoint3d(), c.getPoint3d(), "3D coordinates must match");
        Assertions.assertNotSame(a.getPoint2d(), c.getPoint2d(), "2D point must be a new instance");
        Assertions.assertNotSame(a.getPoint3d(), c.getPoint3d(), "3D point must be a new instance");
    }

    @Test
    void copyPreservesStereochemistryOnCopiedAtoms() throws Exception {
        SmilesParser parser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer original = parser.parseSmiles("C[C@H](N)O"); // tetrahedral centre

        IAtomContainer copy = Descriptor.copyMolecule(original);

        int originalStereoCount = 0;
        for (IStereoElement<?, ?> ignored : original.stereoElements()) {
            originalStereoCount++;
        }
        int copyStereoCount = 0;
        for (IStereoElement<?, ?> element : copy.stereoElements()) {
            copyStereoCount++;
            // Every atom referenced by the stereo element must live in the copy, not the original
            Assertions.assertTrue(copy.contains((IAtom) element.getFocus()),
                    "Stereo focus atom must belong to the copy");
        }
        Assertions.assertTrue(originalStereoCount > 0, "Test input must contain a stereo element");
        Assertions.assertEquals(originalStereoCount, copyStereoCount, "Stereo element count must be preserved");
    }

    /**
     * Test method for getDescriptorAndComponentIndex and getDescriptorAndComponentInfo methods.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_getDescriptorAndComponentIndexAndInfo() throws Exception {
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.BCUT, Descriptor.WIENER_NUMBER};
        // MOLECULAR_WEIGHT (1 component), BCUT (6 components), WIENER_NUMBER (2 components)

        int[] resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 0);
        Assertions.assertArrayEquals(new int[]{0, 0}, resultIndex);
        String[] resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 0);
        Assertions.assertArrayEquals(new String[]{Descriptor.MOLECULAR_WEIGHT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 1);
        Assertions.assertArrayEquals(new int[]{1, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 1);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 6);
        Assertions.assertArrayEquals(new int[]{1, 5}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 6);
        Assertions.assertArrayEquals(new String[]{Descriptor.BCUT.getName(), "5"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 7);
        Assertions.assertArrayEquals(new int[]{2, 0}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 7);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "0"}, resultInfo);

        resultIndex = Descriptor.getDescriptorAndComponentIndex(descriptors, 8);
        Assertions.assertArrayEquals(new int[]{2, 1}, resultIndex);
        resultInfo = Descriptor.getDescriptorAndComponentInfo(descriptors, 8);
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "1"}, resultInfo);

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(descriptors, -1);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(descriptors, -1);
        });

        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentIndex(descriptors, 9);
        });
        Assertions.assertThrows(IllegalArgumentException.class, () -> {
            Descriptor.getDescriptorAndComponentInfo(descriptors, 9);
        });
    }

    /**
     * Test method for setFingerprintPoolSize.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_setFingerprintPoolSize() throws Exception {
        int originalPoolSize = Descriptor.getFingerprintPoolSize();
        try {
            Assertions.assertTrue(Descriptor.setFingerprintPoolSize(6));
            Assertions.assertEquals(6, Descriptor.getFingerprintPoolSize());

            Assertions.assertThrows(IllegalArgumentException.class, () -> {
                Descriptor.setFingerprintPoolSize(0);
            });
            Assertions.assertThrows(IllegalArgumentException.class, () -> {
                Descriptor.setFingerprintPoolSize(-1);
            });
        } finally {
            Descriptor.setFingerprintPoolSize(originalPoolSize);
        }
    }

    /**
     * Test method for setCircularFingerprintSize().
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_setCircularFingerprintSize() throws Exception {
        int originalSize = Descriptor.getCircularFingerprintSize();
        try {
            Assertions.assertTrue(Descriptor.setCircularFingerprintSize(2048));
            Assertions.assertEquals(2048, Descriptor.getCircularFingerprintSize());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6.getDescriptorComponentNumber());
            Assertions.assertEquals(2048, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6.getDescriptorComponentNumber());

            Assertions.assertThrows(IllegalArgumentException.class, () -> {
                Descriptor.setCircularFingerprintSize(0);
            });
            Assertions.assertThrows(IllegalArgumentException.class, () -> {
                Descriptor.setCircularFingerprintSize(-1);
            });
        } finally {
            Descriptor.setCircularFingerprintSize(originalSize);
        }
    }

    /**
     * Sets calculated descriptor components in vectors (rows) of a matrix (that corresponds to atomContainerArray)
     * beginning with startIndex by (optional) parallelization of molecules. If parallel computation is used, the atom
     * containers (molecules) are distributed onto parallel thread, one for each molecule, and they all access shared
     * descriptor instances.
     * Note: Uses a new descriptor instance for EVERY descriptor calculation which slows down the calculation.
     * Note: For fingerprints a blocked queue is used, because fingerprinter instances are not threadsafe.
     * The pool size can be changed via setFingerprintPoolSize(int)
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainerArray Array of molecules. Note: atomContainerArray[i] corresponds to matrix[i] data
     *                             vector. (IS NOT CHANGED)
     * @param matrix Matrix of component vectors of molecules. Note: Data vector matrix[i] corresponds to molecule
     *                atomContainerArray[i], i.e. the molecules define the rows of the matrix. (MAY BE CHANGED)
     * @param startIndex Start index in a vector to be filled with calculated components of descriptors, i.e. matrix
     *                    column to start filling with descriptors
     * @param isParallelCalculation True: Calculations are parallelized, false: Calculations are sequential
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     *                      IMPORTANT: For parallel calculations (isParallelCalculation=true), this must be thread-safe.
     *                      Use Collections.synchronizedList() to avoid race conditions.
     * @return True: Operation was successful, no NaN values generated; false: Operation failed, i.e. at least one component in a descriptor
     *         calculation is NaN or a global exception occurred
     * @throws IllegalArgumentException if the matrix dimensions are invalid
     */
    public static boolean setDescriptorsForMoleculesByMoleculeParallelizationNew(
            Descriptor[] descriptors,
            IAtomContainer[] atomContainerArray,
            float[][] matrix,
            int startIndex,
            boolean isParallelCalculation,
            List<int[]> nanPositionsList
    ) throws IllegalArgumentException {
        // Checks
        final String methodName = "setDescriptorsForMoleculesByMoleculeParallelizationNew";
        if (!Descriptor.validateDescriptors(descriptors, methodName)) {
            Descriptor.LOGGER.warn(methodName + " : Given descriptor array is empty, calculation aborted.");
            return true;
        }
        if (!Descriptor.validateAtomContainerArray(atomContainerArray, methodName)) {
            Descriptor.LOGGER.warn(methodName + " : Given atom container array is empty, calculation aborted.");
            return true;
        }
        for (Descriptor descriptor : descriptors) {
            if (descriptor.requires3DCoordinates()) {
                if (!GeometryUtil.has3DCoordinates(atomContainerArray[0])) {
                    Descriptor.LOGGER.warn(methodName + " : Descriptors require 3D coordinates, but the given molecule does not have 3D coordinates. Calculation aborted.");
                    return true;
                }
                break;
            }
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
                                boolean success = DescriptorTest.setDescriptorsForSingleMoleculeNew(
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
                            } catch (Exception exception) {
                                hasNaN.set(true);
                                Descriptor.LOGGER.warn(
                                        String.format("DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew: One descriptor calculation caused an exception, molecule index: %d.", i),
                                        exception
                                );
                            }
                        }
                );
            } else {
                for (int i = 0; i < atomContainerArray.length; i++) {
                    try {
                        if (!DescriptorTest.setDescriptorsForSingleMoleculeNew(descriptors, atomContainerArray[i], matrix[i], startIndices, i, nanPositionsList)) {
                            hasNaN.set(true);
                        }
                    } catch (Exception exception) {
                        hasNaN.set(true);
                        Descriptor.LOGGER.warn(
                                String.format("DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew: Exception for molecule index: %d", i),
                                exception
                        );
                    }
                }
            }
        } catch (Exception exception) {
            Descriptor.LOGGER.warn(
                    "DescriptorTest.setDescriptorsForMoleculesByMoleculeParallelizationNew: Global exception occurred in descriptor calculation: ",
                    exception
            );
            return false;
        }
        return !hasNaN.get();
    }



    /**
     * Sets calculated descriptor components in vector (that corresponds to atomContainer, a row in the data matrix)
     * at aStartIndices.
     * Note: Uses a new descriptor instance for EVERY descriptor calculation.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptors Array of descriptors to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Component vector of molecule (MAY BE CHANGED)
     * @param aStartIndices Start indices in vector to be filled with calculated components of descriptor (each descriptor
     *                      in descriptors has its dedicated start index here)
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     * @throws CloneNotSupportedException if copying the molecule fails
     */
    private static boolean setDescriptorsForSingleMoleculeNew(
            Descriptor[] descriptors,
            IAtomContainer atomContainer,
            float[] vector,
            int[] aStartIndices,
            int moleculeIndex,
            List<int[]> nanPositionsList
    ) throws CloneNotSupportedException {
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
            if (!DescriptorTest.setDescriptorNew(descriptors[i], moleculeToUse, vector, aStartIndices[i], moleculeIndex, nanPositionsList)) {
                isSuccessful = false;
            }
        }
        return isSuccessful;
    }

    /**
     * Sets component values of descriptor for atomContainer in vector beginning with startIndex.
     * Note: This method instantiates a new CDK descriptor instance for EVERY calculation and is therefore thread-safe in
     * concurrent computing.
     * Note: Checks are NOT performed here. All necessary checks have already been made in public methods above.
     *
     * @param descriptor Descriptor to be calculated (IS NOT CHANGED)
     * @param atomContainer Molecule (IS NOT CHANGED)
     * @param vector Vector of molecule (row in data matrix) to be filled with calculated components of descriptors (MAY BE CHANGED)
     * @param startIndex Start index in vector to be filled with calculated components of descriptors
     * @param moleculeIndex Index of the current molecule being processed
     * @param nanPositionsList List to track NaN positions as [moleculeIndex, componentIndex] pairs (MAY BE CHANGED).
     * @return True: Operation was successful, no NaN values were generated; false: Operation failed, i.e. at least one component in a
     * descriptor calculation result is NaN
     */
    private static boolean setDescriptorNew(
            Descriptor descriptor,
            IAtomContainer atomContainer,
            float[] vector,
            int startIndex,
            int moleculeIndex,
            List<int[]> nanPositionsList
    ) {
        try {
            switch (descriptor) {
                //note: the initializations here should be the same as in the descriptorToCdkObjectMap!
                case MOLECULAR_WEIGHT:
                    vector[startIndex] = (float) ((DoubleResult) (new WeightDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case WIENER_NUMBER:
                    DoubleArrayResult result = (DoubleArrayResult) (new WienerNumbersDescriptor()).calculate(atomContainer).getValue();
                    vector[startIndex] = (float) result.get(0); //Wiener path number
                    vector[startIndex + 1] = (float) result.get(1); //Wiener polarity number
                    break;
                case ATOM_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new AtomCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_HEAVY:
                    AtomCountDescriptor atomCountHeavyDesc = new AtomCountDescriptor();
                    //set parameter to count heavy atoms
                    atomCountHeavyDesc.setParameters(new String[]{"#"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountHeavyDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_C:
                    AtomCountDescriptor atomCountCDesc = new AtomCountDescriptor();
                    //set parameter to count carbon atoms
                    atomCountCDesc.setParameters(new String[]{"C"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountCDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_H:
                    AtomCountDescriptor atomCountHDesc = new AtomCountDescriptor();
                    //set parameter to count hydrogen atoms
                    atomCountHDesc.setParameters(new String[]{"H"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountHDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_N:
                    AtomCountDescriptor atomCountNDesc = new AtomCountDescriptor();
                    //set parameter to count nitrogen atoms
                    atomCountNDesc.setParameters(new String[]{"N"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountNDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_O:
                    AtomCountDescriptor atomCountODesc = new AtomCountDescriptor();
                    //set parameter to count oxygen atoms
                    atomCountODesc.setParameters(new String[]{"O"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountODesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_S:
                    AtomCountDescriptor atomCountSDesc = new AtomCountDescriptor();
                    //set parameter to count sulfur atoms
                    atomCountSDesc.setParameters(new String[]{"S"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountSDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_P:
                    AtomCountDescriptor atomCountPDesc = new AtomCountDescriptor();
                    //set parameter to count phosphorus atoms
                    atomCountPDesc.setParameters(new Object[]{"P"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountPDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_F:
                    AtomCountDescriptor atomCountFDesc = new AtomCountDescriptor();
                    //set parameter to count fluorine atoms
                    atomCountFDesc.setParameters(new String[]{"F"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountFDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_BR:
                    AtomCountDescriptor atomCountBrDesc = new AtomCountDescriptor();
                    //set parameter to count bromine atoms
                    atomCountBrDesc.setParameters(new String[]{"Br"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountBrDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_CL:
                    AtomCountDescriptor atomCountClDesc = new AtomCountDescriptor();
                    //set parameter to count chlorine atoms
                    atomCountClDesc.setParameters(new String[]{"Cl"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountClDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ATOM_COUNT_I:
                    AtomCountDescriptor atomCountIDesc = new AtomCountDescriptor();
                    //set parameter to count iodine atoms
                    atomCountIDesc.setParameters(new String[]{"I"});
                    vector[startIndex] = (float) ((IntegerResult) atomCountIDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case H_BOND_ACCEPTOR_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new HBondAcceptorCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case H_BOND_DONOR_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new HBondDonorCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case TPSA:
                    vector[startIndex] = (float) ((DoubleResult) (new TPSADescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_CHAIN:
                    vector[startIndex] = (float) ((IntegerResult) (new LargestChainDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case LONGEST_ALIPHATIC_CHAIN:
                    vector[startIndex] = (float) ((IntegerResult) (new LongestAliphaticChainDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case MANNHOLD_LOGP:
                    vector[startIndex] = (float) ((DoubleResult) (new MannholdLogPDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case BCUT:
                    BCUTDescriptor bcutDescriptor = new BCUTDescriptor();
                    // Change parameters so that we can use our own setAromaticity method default: checkAromaticity = true
                    bcutDescriptor.setParameters(new Object[] {1, 1, false}); // nhigh = 1, nlow = 1, checkAromaticity = false
                    DoubleArrayResult bcutResult = (DoubleArrayResult) bcutDescriptor.calculate(atomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        vector[startIndex + i] = (float) bcutResult.get(i);
                    }
                    break;
                case BOND_COUNT_ALL:
                    vector[startIndex] = (float) ((IntegerResult) (new BondCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_SINGLE:
                    BondCountDescriptor bondCountSingleDesc = new BondCountDescriptor();
                    //set parameter to count single bonds
                    bondCountSingleDesc.setParameters(new String[]{"s"});
                    vector[startIndex] = (float) ((IntegerResult) bondCountSingleDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_DOUBLE:
                    BondCountDescriptor bondCountDoubleDesc = new BondCountDescriptor();
                    //set parameter to count double bonds
                    bondCountDoubleDesc.setParameters(new String[]{"d"});
                    vector[startIndex] = (float) ((IntegerResult) bondCountDoubleDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case BOND_COUNT_TRIPLE:
                    BondCountDescriptor bondCountTripleDesc = new BondCountDescriptor();
                    //set parameter to count triple bonds
                    bondCountTripleDesc.setParameters(new String[]{"t"});
                    vector[startIndex] = (float) ((IntegerResult) bondCountTripleDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case B_POL:
                    vector[startIndex] = (float) ((DoubleResult) (new BPolDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case RULE_OF_FIVE:
                    vector[startIndex] = (float) ((IntegerResult) (new RuleOfFiveDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case AROMATIC_ATOMS_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new AromaticAtomsCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case AROMATIC_BONDS_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new AromaticBondsCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case ROTATABLE_BONDS_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new RotatableBondsCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case FMF:
                    vector[startIndex] = (float) ((DoubleResult) (new FMFDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case FRACTIONAL_CSP3:
                    vector[startIndex] = (float) ((DoubleResult) (new FractionalCSP3Descriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case HYBRIDIZATION_RATIO:
                    vector[startIndex] = (float) ((DoubleResult) (new HybridizationRatioDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case KAPPA_SHAPE_INDICES:
                    DoubleArrayResult kappaResult = (DoubleArrayResult) (new KappaShapeIndicesDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 3; i++) {
                        vector[startIndex + i] = (float) kappaResult.get(i);
                    }
                    break;
                case PETITJEAN_NUMBER:
                    vector[startIndex] = (float) ((DoubleResult) (new PetitjeanNumberDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case SPIRO_ATOM_COUNT:
                    vector[startIndex] = (float) ((IntegerResult) (new SpiroAtomCountDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case V_ADJ_MAT:
                    vector[startIndex] = (float) ((DoubleResult) (new VAdjMaDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case WEIGHTED_PATH:
                    DoubleArrayResult weightedPathResultNew = (DoubleArrayResult) (new WeightedPathDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        vector[startIndex + i] = (float) weightedPathResultNew.get(i);
                    }
                    break;
                case ZAGREB_INDEX:
                    vector[startIndex] = (float) ((DoubleResult) (new ZagrebIndexDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case CARBON_TYPES:
                    IntegerArrayResult carbonTypesResultNew = (IntegerArrayResult) (new CarbonTypesDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 9; i++) {
                        vector[startIndex + i] = carbonTypesResultNew.get(i);
                    }
                    break;
                case A_LOG_P:
                    DoubleArrayResult aLogPResult = (DoubleArrayResult) (new ALOGPDescriptor()).calculate(atomContainer).getValue();
                    vector[startIndex] = (float) aLogPResult.get(0);// ALogP
                    vector[startIndex + 1] = (float) aLogPResult.get(1);  // ALogP squared
                    vector[startIndex + 2] = (float) aLogPResult.get(2);  // Molar Refractivity
                    break;
                case X_LOG_P:
                    vector[startIndex] = (float) ((DoubleResult) (new XLogPDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case JP_LOG_P:
                    vector[startIndex] = (float) ((DoubleResult) (new JPlogPDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case A_POL:
                    vector[startIndex] = (float) ((DoubleResult) (new APolDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case AUTOCORRELATION_CHARGE:
                    DoubleArrayResult autocorrelationChargeResult = (DoubleArrayResult) (new AutocorrelationDescriptorCharge()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        vector[startIndex + i] = (float) autocorrelationChargeResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_MASS:
                    DoubleArrayResult autocorrelationMassResult = (DoubleArrayResult) (new AutocorrelationDescriptorMass()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        vector[startIndex + i] = (float) autocorrelationMassResult.get(i);
                    }
                    break;
                case AUTOCORRELATION_POLARIZABILITY:
                    DoubleArrayResult autocorrelationPolarizabilityResult = (DoubleArrayResult) (new AutocorrelationDescriptorPolarizability()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 5; i++) {
                        vector[startIndex + i] = (float) autocorrelationPolarizabilityResult.get(i);
                    }
                    break;
                case FRAGMENT_COMPLEXITY:
                    vector[startIndex] = (float) ((DoubleResult) (new FragmentComplexityDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case CHI_CHAIN:
                    DoubleArrayResult chiChainResult = (DoubleArrayResult) new ChiChainDescriptor().calculate(atomContainer).getValue();
                    for (int i = 0; i < 10; i++) {
                        vector[startIndex + i] = (float) chiChainResult.get(i);
                    }
                    break;
                case CHI_CLUSTER:
                    DoubleArrayResult chiClusterResult = (DoubleArrayResult) (new ChiClusterDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 8; i++) {
                        vector[startIndex + i] = (float) chiClusterResult.get(i);
                    }
                    break;
                case CHI_PATH_CLUSTER:
                    DoubleArrayResult chiPathClusterResult = (DoubleArrayResult) new ChiPathClusterDescriptor().calculate(atomContainer).getValue();
                    for (int i = 0; i < 6; i++) {
                        vector[startIndex + i] = (float) chiPathClusterResult.get(i);
                    }
                    break;
                case CHI_PATH:
                    DoubleArrayResult chiPathResult = (DoubleArrayResult) (new ChiPathDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 16; i++) {
                        vector[startIndex + i] = (float) chiPathResult.get(i);
                    }
                    break;
                case FRACTIONAL_PSA:
                    vector[startIndex] = (float) ((DoubleResult) (new FractionalPSADescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case LARGEST_PI_SYSTEM:
                    LargestPiSystemDescriptor largestPiSystemDescriptor = new LargestPiSystemDescriptor();
                    // Change parameters so that we can use our own setAromaticity method
                    largestPiSystemDescriptor.setParameters(new Boolean[] {false}); // checkAromaticity = false
                    vector[startIndex] = (float) ((IntegerResult) largestPiSystemDescriptor.calculate(atomContainer).getValue()).intValue();
                    break;
                case SMALL_RING:
                    IntegerArrayResult smallRingResult = (IntegerArrayResult) (new SmallRingDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 11; i++) {
                        vector[startIndex + i] = (float) smallRingResult.get(i);
                    }
                    break;
                case BASIC_GROUP_COUNT:
                    BasicGroupCountDescriptor basicGroupCountDesc = new BasicGroupCountDescriptor();
                    basicGroupCountDesc.initialise(SilentChemObjectBuilder.getInstance());
                    vector[startIndex] = (float) ((IntegerResult) basicGroupCountDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case ACIDIC_GROUP_COUNT:
                    AcidicGroupCountDescriptor acidicGroupCountDesc = new AcidicGroupCountDescriptor();
                    acidicGroupCountDesc.initialise(SilentChemObjectBuilder.getInstance());
                    vector[startIndex] = (float) ((IntegerResult) acidicGroupCountDesc.calculate(atomContainer).getValue()).intValue();
                    break;
                case AMINO_ACID_COUNT:
                    IntegerArrayResult aminoAcidCountResult = (IntegerArrayResult) (new AminoAcidCountDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 20; i++) {
                        vector[startIndex + i] = (float) aminoAcidCountResult.get(i);
                    }
                    break;
                case KIER_HALL_SMARTS:
                    IntegerArrayResult kierHallSmartsResult = (IntegerArrayResult) (new KierHallSmartsDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 79; i++) {
                        vector[startIndex + i] = (float) kierHallSmartsResult.get(i);
                    }
                    break;
                case ECCENTRIC_CONNECTIVITY_INDEX:
                    vector[startIndex] = (float) ((IntegerResult) (new EccentricConnectivityIndexDescriptor()).calculate(atomContainer).getValue()).intValue();
                    break;
                case MDE:
                    DoubleArrayResult aDoubleArrayResult = (DoubleArrayResult) (new MDEDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 19; i++) {
                        vector[startIndex + i] = (float) aDoubleArrayResult.get(i);
                    }
                    break;
                case VABC:
                    vector[startIndex] = (float) ((DoubleResult) (new VABCDescriptor()).calculate(atomContainer).getValue()).doubleValue();
                    break;
                case PUBCHEM_FINGERPRINTER:
                    //TODO: why do the fingerprints have their own try-catch block? It does the same as the general one below, right? So, can we remove them?
                    try {
                        IBitFingerprint fingerprint = new PubchemFingerprinter(SilentChemObjectBuilder.getInstance()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.PUBCHEM_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.PUBCHEM_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_0:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP0, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_0:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP0, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_2:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP2, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_2:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP2, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_4:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_4:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP4, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_ECFP_6:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP6, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case CIRCULAR_FINGERPRINTER_FCFP_6:
                    try {
                        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_FCFP6, Descriptor.getCircularFingerprintSize()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.getCircularFingerprintSize(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;
                case MACCS_FINGERPRINTER:
                    try {
                        IBitFingerprint fingerprint = new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()).getBitFingerprint(atomContainer);
                        for (int i = 0; i < Descriptor.MACCS_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = fingerprint.get(i) ? 1.0f : 0.0f;
                        }
                    } catch (Exception exception) {
                        for (int i = 0; i < Descriptor.MACCS_FINGERPRINTER.getDescriptorComponentNumber(); i++) {
                            vector[startIndex + i] = Float.NaN;
                        }
                        Descriptor.LOGGER.warn(exception.toString(), exception);
                    }
                    break;

                case CPSA:
                    DoubleArrayResult cpsaResult = (DoubleArrayResult) (new CPSADescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 29; i++) {
                        vector[startIndex + i] = (float) cpsaResult.get(i);
                    }
                    break;
                case GRAVITATIONAL_INDEX:
                    DoubleArrayResult gravResult = (DoubleArrayResult) (new GravitationalIndexDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 9; i++) {
                        vector[startIndex + i] = (float) gravResult.get(i);
                    }
                    break;
                case MOMENT_OF_INERTIA:
                    DoubleArrayResult momiResult = (DoubleArrayResult) (new MomentOfInertiaDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 7; i++) {
                        vector[startIndex + i] = (float) momiResult.get(i);
                    }
                    break;
                case WHIM:
                    DoubleArrayResult whimResult = (DoubleArrayResult) (new WHIMDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 17; i++) {
                        vector[startIndex + i] = (float) whimResult.get(i);
                    }
                    break;
                case LENGTH_OVER_BREADTH:
                    DoubleArrayResult lobResult = (DoubleArrayResult) (new LengthOverBreadthDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 2; i++) {
                        vector[startIndex + i] = (float) lobResult.get(i);
                    }
                    break;
                case PETITJEAN_SHAPE_INDEX:
                    DoubleArrayResult pShapeResult = (DoubleArrayResult) (new PetitjeanShapeIndexDescriptor()).calculate(atomContainer).getValue();
                    for (int i = 0; i < 2; i++) {
                        vector[startIndex + i] = (float) pShapeResult.get(i);
                    }
                    break;

                // Add new descriptor information here!
                default:
                    throw new UnsupportedOperationException(descriptor + ": This descriptor does not have a routine yet!");
            }
            // Check for NaN values in the calculated result and track them
            int numComponents = descriptor.getDescriptorComponentNumber();
            return !Descriptor.checkAndTrackNaNValues(vector, startIndex, numComponents, moleculeIndex, nanPositionsList);
        } catch (Exception exception) {
            int numComponents = descriptor.getDescriptorComponentNumber();
            for (int i = 0; i < numComponents; i++) {
                vector[startIndex + i] = Float.NaN;
                // Track NaN position if nanPositionsList is provided
                if (nanPositionsList != null) {
                    nanPositionsList.add(new int[]{moleculeIndex, startIndex + i});
                }
            }
            Descriptor.LOGGER.warn(
                    String.format("DescriptorTest.setDescriptorNew: An exception occurred while calculating descriptor %s for molecule index %d.",
                            descriptor,
                            moleculeIndex),
                    exception
            );
            return false;
        }
    }

    @Test
    void test_GetAllFingerprints_Excludes3D() {
        Descriptor[] fingerprints = Descriptor.getAllFingerprints();
        Assertions.assertTrue(fingerprints.length > 0);
        for (Descriptor fp : fingerprints) {
            Assertions.assertTrue(fp.isFingerprint());
            Assertions.assertFalse(fp.requires3DCoordinates());
        }
    }

    @Test
    void test_GetSpecifiedDescriptors_3DInclusion() {
        Descriptor[] without3D = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        for (Descriptor d : without3D) {
            Assertions.assertFalse(d.requires3DCoordinates());
        }

        Descriptor[] with3D = Descriptor.getSpecifiedDescriptors(true, true, true, true);
        Assertions.assertEquals(Descriptor.values().length, with3D.length);

        Descriptor[] defaultSelection = Descriptor.getSpecifiedDescriptors(false, false, false, false);
        for (Descriptor d : defaultSelection) {
            Assertions.assertTrue(d.isFast());
            Assertions.assertTrue(d.isSafe());
            Assertions.assertFalse(d.isFingerprint());
            Assertions.assertFalse(d.requires3DCoordinates());
        }
    }

    @Test
    void test_Requires3DCoordinates_Getter() {
        for (Descriptor descriptor : Descriptor.values()) {
            boolean is3D = descriptor == Descriptor.CPSA
                    || descriptor == Descriptor.GRAVITATIONAL_INDEX
                    || descriptor == Descriptor.MOMENT_OF_INERTIA
                    || descriptor == Descriptor.WHIM
                    || descriptor == Descriptor.LENGTH_OVER_BREADTH
                    || descriptor == Descriptor.PETITJEAN_SHAPE_INDEX;
            Assertions.assertEquals(is3D, descriptor.requires3DCoordinates(), "Mismatch for " + descriptor);
        }
    }
}
