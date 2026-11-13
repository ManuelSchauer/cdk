/*
 * Copyright (C) 2026 Manuel Schauer, Jonas Schaub, Christoph Steinbeck, and Achim Zielesny
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */

package org.openscience.cdk.qsar;


import org.junit.jupiter.api.Assertions;
import org.junit.jupiter.api.Tag;
import org.junit.jupiter.api.Test;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.fingerprint.CircularFingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.MACCSFingerprinter;
import org.openscience.cdk.fingerprint.PubchemFingerprinter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorCharge;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorMass;
import org.openscience.cdk.qsar.descriptors.molecular.AutocorrelationDescriptorPolarizability;
import org.openscience.cdk.qsar.descriptors.molecular.PetitjeanShapeIndexDescriptor;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.config.Isotopes;
import org.openscience.cdk.io.HINReader;
import org.openscience.cdk.io.ISimpleChemObjectReader;
import org.openscience.cdk.io.MDLV2000Reader;
import org.openscience.cdk.qsar.descriptors.molecular.CPSADescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.GravitationalIndexDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.LengthOverBreadthDescriptor;
import org.openscience.cdk.qsar.descriptors.molecular.MomentOfInertiaDescriptor;
import org.openscience.cdk.qsar.result.DoubleArrayResult;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import java.io.InputStream;

import javax.vecmath.Point2d;
import javax.vecmath.Point3d;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

/**
 * Test class for DescriptorCalculator class. This class first tests all the descriptors included in the {@link Descriptor} class
 * individually, i.e. whether they produce the expected results for some example molecules (in most cases taken from the
 * respective CDK descriptor test class) in different parallelization settings.
 * <p>
 * Note on the scope of these tests: they check that the engine calculates the CDK descriptors and fingerprints
 * correctly, i.e. that its results are equal to reference values from the CDK descriptor tests or, where no reference
 * values exist (e.g. autocorrelation descriptors, most circular fingerprints), to the results of the CDK classes
 * calculated directly, without the engine. They do NOT test whether the descriptors themselves are scientifically
 * correct; this is the responsibility of the respective CDK descriptor and fingerprint tests.
 * <p>
 * Note: For adding tests of a new descriptor goto "Add new descriptor tests here!"
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorCalculatorTest {
    /**
     * Calculator with the default configuration, shared by all tests (a calculator can be used by any number of
     * threads at the same time).
     */
    private static final DescriptorCalculator CALCULATOR = new DescriptorCalculator();

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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            DescriptorCalculator.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            DescriptorCalculator.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            DescriptorCalculator.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        // Apply aromaticity to the molecules (LargestChainDescriptor does not depend on aromaticity models)
        DescriptorCalculator.setAromaticity(molecule1, Aromaticity.Model.Daylight);
        DescriptorCalculator.setAromaticity(molecule2, Aromaticity.Model.Daylight);

        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule1, model);
            // Apply aromaticity with the current model to the second molecule
            DescriptorCalculator.setAromaticity(molecule2, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule1, molecule2};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
     * <p>
     * There are no validated reference values for this descriptor. Therefore, the result of the engine is compared with
     * the result of the CDK class
     * {@link AutocorrelationDescriptorCharge} calculated directly, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_CHARGE() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        // reference: the CDK descriptor calculated directly, without the engine (it works on a copy of the molecule)
        float[] expected = DescriptorCalculatorTest.toFloatArray((DoubleArrayResult) new AutocorrelationDescriptorCharge().calculate(molecule).getValue());
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_CHARGE};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0], "engine must equal the direct CDK calculation");
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_MASS.
     * <p>
     * There are no validated reference values for this descriptor. Therefore, the result of the engine is compared with
     * the result of the CDK class
     * {@link AutocorrelationDescriptorMass} calculated directly, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_MASS() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        // reference: the CDK descriptor calculated directly, without the engine (it works on a copy of the molecule)
        float[] expected = DescriptorCalculatorTest.toFloatArray((DoubleArrayResult) new AutocorrelationDescriptorMass().calculate(molecule).getValue());
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_MASS};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0], "engine must equal the direct CDK calculation");
        Assertions.assertNotEquals(0f, matrix[0][0]);
        Assertions.assertNotEquals(0f, matrix[0][1]);
        Assertions.assertNotEquals(0f, matrix[0][2]);
        Assertions.assertNotEquals(0f, matrix[0][3]);
        Assertions.assertNotEquals(0f, matrix[0][4]);
    }

    /**
     * Test method for descriptor AUTOCORRELATION_POLARIZABILITY.
     * <p>
     * There are no validated reference values for this descriptor. Therefore, the result of the engine is compared with
     * the result of the CDK class
     * {@link AutocorrelationDescriptorPolarizability} calculated directly, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_AUTOCORRELATION_POLARIZABILITY() throws Exception {
        String smiles = "Clc1ccccc1"; // Chlorobenzene CID: 7964
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles(smiles);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        // reference: the CDK descriptor calculated directly, without the engine (it works on a copy of the molecule)
        float[] expected = DescriptorCalculatorTest.toFloatArray((DoubleArrayResult) new AutocorrelationDescriptorPolarizability().calculate(molecule).getValue());
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.AUTOCORRELATION_POLARIZABILITY};
        boolean isParallelCalculation = false;

        Assertions.assertEquals(5, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        Assertions.assertArrayEquals(expected, matrix[0], "engine must equal the direct CDK calculation");
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
            DescriptorCalculator.setAromaticity(molecule, model);

            IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};

            float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
            Assertions.assertTrue(
                    CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrix,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );
        // Note: Glycine (index 8) is matched for all alpha-amino acids because the glycine substructure/SMARTS
        // matches the common amino acid backbone core (inherent CDK AminoAcidCountDescriptor behavior).
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
        DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        int startIndex = 0;
        Descriptor[] descriptors = new Descriptor[]{Descriptor.VABC};
        boolean isParallelCalculation = false;
        double epsilon = 0.01; // tolerance range
        float expected = 292.23f;

        Assertions.assertEquals(1, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[moleculesArray.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
        DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
     * <p>
     * Reference values: the CDK test class MACCSFingerprinterTest has no complete reference fingerprint, only single
     * reference bits, which are tested here: 1,2-diphenylethane sets bit 124 but not bit 165, ethylbenzene does not set
     * bit 124, and "CCC.CCC" (two fragments) sets bit 165. In addition, all 166 bits are compared with the fingerprint
     * calculated directly with the CDK class {@link MACCSFingerprinter}, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_MACCS_FINGERPRINTER() throws Exception {
        String[] smilesArray = new String[]{
                "c1ccccc1CCc1ccccc1", // 1,2-Diphenylethane CID: 7647
                "c1ccccc1CC", // Ethylbenzene CID: 7500
                "CCC.CCC" // two propane fragments
        };
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer[] moleculesArray = new IAtomContainer[smilesArray.length];
        float[][] expected = new float[smilesArray.length][];
        for (int i = 0; i < smilesArray.length; i++) {
            moleculesArray[i] = smilesParser.parseSmiles(smilesArray[i]);
            DescriptorCalculator.setAromaticity(moleculesArray[i], Aromaticity.Model.Daylight);
            // reference: the fingerprint calculated directly with the CDK class, without the engine
            expected[i] = DescriptorCalculatorTest.toFloatArray(
                    new MACCSFingerprinter(SilentChemObjectBuilder.getInstance()).getBitFingerprint(moleculesArray[i]), 166);
        }
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MACCS_FINGERPRINTER};
        Assertions.assertEquals(166, Descriptor.getNumberOfComponents(descriptors));

        float[][] matrix = new float[moleculesArray.length][166];
        List<int[]> nanPositions = new ArrayList<>();
        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, moleculesArray, matrix, 0, false, nanPositions));
        // reference bits taken from MACCSFingerprinterTest
        Assertions.assertEquals(1f, matrix[0][124]);
        Assertions.assertEquals(0f, matrix[0][165]);
        Assertions.assertEquals(0f, matrix[1][124]);
        Assertions.assertEquals(1f, matrix[2][165]);
        // complete comparison with the direct CDK calculation
        for (int i = 0; i < moleculesArray.length; i++) {
            Assertions.assertArrayEquals(expected[i], matrix[i], smilesArray[i]);
        }
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_ECFP descriptors of different diameters.
     * <p>
     * Reference values: for ECFP6 with 1024 bits, the CDK test class CircularFingerprinterTest contains a complete
     * reference fingerprint of propanamide ("CCC(=O)N", bits 19, 152, 293, 340, 439, 480, 507, 726, 762, 947, 993),
     * which is tested here. For ECFP0, ECFP2 and ECFP4 there are no reference fingerprints; therefore, all bits of all
     * four fingerprints are additionally compared with the fingerprints calculated directly with the CDK class
     * {@link CircularFingerprinter}, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_ECFP() throws Exception {
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6
        };
        int[] cdkClasses = new int[]{
                CircularFingerprinter.CLASS_ECFP0,
                CircularFingerprinter.CLASS_ECFP2,
                CircularFingerprinter.CLASS_ECFP4,
                CircularFingerprinter.CLASS_ECFP6
        };
        // row 1 is propanamide, ECFP6 is the fourth fingerprint in the row
        final int[] referenceBitsEcfp6 = new int[]{19, 152, 293, 340, 439, 480, 507, 726, 762, 947, 993};
        int size = CALCULATOR.getCircularFingerprintSize();
        Assertions.assertEquals(1024, size, "reference fingerprint requires the default size of 1024 bits");
        float[][] matrix = DescriptorCalculatorTest.calculateCircularFingerprintsAndCompareWithCdk(descriptors, cdkClasses);
        BitSet expectedBits = new BitSet();
        for (int bit : referenceBitsEcfp6) {
            expectedBits.set(bit);
        }
        for (int i = 0; i < size; i++) {
            Assertions.assertEquals(expectedBits.get(i) ? 1f : 0f, matrix[1][3 * size + i], "ECFP6 bit " + i);
        }
    }

    /**
     * Test method for all CIRCULAR_FINGERPRINTER_FCFP descriptors of different diameters.
     * <p>
     * There are no reference fingerprints for FCFP in the CDK test classes (only raw FCFP6 features in the
     * CircularFingerprinterTest validation data, not folded bit fingerprints). Therefore, all bits of all four
     * fingerprints are compared with the fingerprints calculated directly with the CDK class
     * {@link CircularFingerprinter}, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CIRCULAR_FINGERPRINTER_FCFP() throws Exception {
        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6
        };
        int[] cdkClasses = new int[]{
                CircularFingerprinter.CLASS_FCFP0,
                CircularFingerprinter.CLASS_FCFP2,
                CircularFingerprinter.CLASS_FCFP4,
                CircularFingerprinter.CLASS_FCFP6
        };
        DescriptorCalculatorTest.calculateCircularFingerprintsAndCompareWithCdk(descriptors, cdkClasses);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[] singleResult = CALCULATOR.calculateDescriptor(Descriptor.CPSA, molecule);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[] singleResult = CALCULATOR.calculateDescriptor(Descriptor.GRAVITATIONAL_INDEX, molecule);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[] singleResult = CALCULATOR.calculateDescriptor(Descriptor.MOMENT_OF_INERTIA, molecule);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[] singleResult = CALCULATOR.calculateDescriptor(Descriptor.LENGTH_OVER_BREADTH, molecule);
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
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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

        float[] singleResult = CALCULATOR.calculateDescriptor(Descriptor.PETITJEAN_SHAPE_INDEX, molecule);
        Assertions.assertEquals(0.5f, singleResult[0], epsilon);
        Assertions.assertEquals(0.606477f, singleResult[1], epsilon);
    }

    // Add new descriptor tests here!

    /**
     * Calculates the given circular fingerprints for 1,2-diphenylethane (with Daylight aromaticity) and propanamide
     * (parsed only, as in the CDK test class CircularFingerprinterTest) with the engine and asserts that all bits are equal to the fingerprints calculated directly
     * with the CDK class {@link CircularFingerprinter}, i.e. without the engine.
     * <p>
     * Note: This test checks that the engine calculates the descriptor correctly (correct CDK class and parameters,
     * correct conversion of the result and correct position in the data matrix), NOT that the descriptor itself is
     * scientifically correct; the latter is the responsibility of the respective CDK descriptor tests.
     *
     * @param descriptors the circular fingerprint descriptors to test
     * @param cdkClasses the corresponding CircularFingerprinter class constants (same order as descriptors)
     * @return the result matrix of the engine (row 0: 1,2-diphenylethane, row 1: propanamide) for additional checks
     *         against reference values
     * @throws Exception if anything goes wrong
     */
    private static float[][] calculateCircularFingerprintsAndCompareWithCdk(Descriptor[] descriptors, int[] cdkClasses) throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer diphenylethane = smilesParser.parseSmiles("c1ccccc1CCc1ccccc1"); // 1,2-Diphenylethane CID: 7647
        DescriptorCalculator.setAromaticity(diphenylethane, Aromaticity.Model.Daylight);
        IAtomContainer propanamide = smilesParser.parseSmiles("CCC(=O)N"); // Propanamide CID: 6578
        IAtomContainer[] moleculesArray = new IAtomContainer[]{diphenylethane, propanamide};
        int size = CALCULATOR.getCircularFingerprintSize();
        Assertions.assertEquals(descriptors.length * size, Descriptor.getNumberOfComponents(descriptors));

        // reference: the fingerprints calculated directly with the CDK class, without the engine
        float[][] expected = new float[moleculesArray.length][descriptors.length * size];
        for (int m = 0; m < moleculesArray.length; m++) {
            for (int d = 0; d < descriptors.length; d++) {
                float[] bits = DescriptorCalculatorTest.toFloatArray(
                        new CircularFingerprinter(cdkClasses[d], size).getBitFingerprint(moleculesArray[m]), size);
                System.arraycopy(bits, 0, expected[m], d * size, size);
            }
        }

        float[][] matrix = new float[moleculesArray.length][descriptors.length * size];
        List<int[]> nanPositions = new ArrayList<>();
        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, moleculesArray, matrix, 0, false, nanPositions));
        Assertions.assertTrue(nanPositions.isEmpty());
        for (int m = 0; m < moleculesArray.length; m++) {
            Assertions.assertArrayEquals(expected[m], matrix[m], "molecule " + m);
        }
        return matrix;
    }

    /**
     * Converts a CDK DoubleArrayResult to the float values the engine writes into the data matrix.
     *
     * @param result the CDK result
     * @return the values cast to float
     */
    private static float[] toFloatArray(DoubleArrayResult result) {
        float[] values = new float[result.length()];
        for (int i = 0; i < values.length; i++) {
            values[i] = (float) result.get(i);
        }
        return values;
    }

    /**
     * Converts a CDK bit fingerprint to the 0/1 float values the engine writes into the data matrix.
     *
     * @param fingerprint the CDK bit fingerprint
     * @param size the number of bits
     * @return 1.0f for set bits, 0.0f otherwise
     */
    private static float[] toFloatArray(IBitFingerprint fingerprint, int size) {
        float[] values = new float[size];
        for (int i = 0; i < size; i++) {
            values[i] = fingerprint.get(i) ? 1.0f : 0.0f;
        }
        return values;
    }

    /**
     * Helper method to generate an array of independent propionic acid molecule instances with Daylight aromaticity configured.
     *
     * @param numberOfMolecules number of molecule instances to create
     * @return array of {@link IAtomContainer} instances
     * @throws Exception if parsing fails
     */
    private static IAtomContainer[] createPropionicAcidMolecules(int numberOfMolecules) throws Exception {
        String smiles = "CCC(=O)O"; // Propionic acid CID: 1032
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer[] moleculesArray = new IAtomContainer[numberOfMolecules];
        for (int i = 0; i < numberOfMolecules; i++) {
            IAtomContainer molecule = smilesParser.parseSmiles(smiles);
            DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
            moleculesArray[i] = molecule;
        }
        return moleculesArray;
    }

    /**
     * Helper method to generate an array of propionic acid SMILES strings.
     *
     * @param numberOfMolecules number of SMILES strings to generate
     * @return array of SMILES strings
     */
    private static String[] createPropionicAcidSmiles(int numberOfMolecules) {
        String smiles = "CCC(=O)O"; // Propionic acid CID: 1032
        String[] moleculeStringsArray = new String[numberOfMolecules];
        Arrays.fill(moleculeStringsArray, smiles);
        return moleculeStringsArray;
    }

    /**
     * Tests parallelization using {@link DescriptorCalculator#setDescriptorsForMoleculesByMoleculeParallelization}.
     * Verifies that distributing individual {@link IAtomContainer} objects across threads in parallel
     * yields identical results to sequential execution across 1000 molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    @Tag("SlowTest")
    void test_Parallelization_MoleculeParallelization() throws Exception {
        int numberOfMolecules = 1000;
        IAtomContainer[] moleculesArray = createPropionicAcidMolecules(numberOfMolecules);
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);

        // Sequential calculation
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        // Parallel calculation
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        // Compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests parallelization using {@link DescriptorCalculator#setDescriptorsForMoleculesByBatchParallelization}.
     * Verifies that distributing batches of {@link IAtomContainer} objects across threads (batch size 100)
     * yields identical results between sequential and parallel execution across 1000 molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    @Tag("SlowTest")
    void test_Parallelization_BatchParallelization() throws Exception {
        int numberOfMolecules = 1000;
        IAtomContainer[] moleculesArray = createPropionicAcidMolecules(numberOfMolecules);
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        int batchSize = 100;

        // Sequential calculation
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositions
                )
        );

        // Parallel calculation
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors,
                        moleculesArray,
                        matrixParallel,
                        startIndex,
                        batchSize,
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        // Compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests parallelization using {@link DescriptorCalculator#setDescriptorsForMoleculeBySmilesStringsBatchParallelization}.
     * Verifies that distributing batches of SMILES strings across threads (batch size 100) with on-the-fly parsing
     * yields identical results between sequential and parallel execution across 1000 molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    @Tag("SlowTest")
    void test_Parallelization_SmilesStringsBatchParallelization() throws Exception {
        int numberOfMolecules = 1000;
        String[] moleculeStringsArray = createPropionicAcidSmiles(numberOfMolecules);
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        int batchSize = 100;

        // Sequential calculation
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixSequential,
                        startIndex,
                        batchSize,
                        null, // use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );

        // Parallel calculation
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixParallel,
                        startIndex,
                        batchSize,
                        null, // use default aromaticity model
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        // Compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests parallelization using {@link DescriptorCalculator#setDescriptorsForMoleculesBySmilesStringParallelization}.
     * Verifies that distributing individual SMILES strings across threads with on-the-fly parsing
     * yields identical results between sequential and parallel execution across 1000 molecules.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    @Tag("SlowTest")
    void test_Parallelization_SmilesStringParallelization() throws Exception {
        int numberOfMolecules = 1000;
        String[] moleculeStringsArray = createPropionicAcidSmiles(numberOfMolecules);
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);

        // Sequential calculation
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixSequential,
                        startIndex,
                        null, // use default aromaticity model
                        isParallelCalculation,
                        nanPositions
                )
        );

        // Parallel calculation
        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                        descriptors,
                        moleculeStringsArray,
                        matrixParallel,
                        startIndex,
                        null, // use default aromaticity model
                        isParallelCalculation,
                        nanPositionsParallel
                )
        );

        // Compare the results
        for (int i = 0; i < numberOfMolecules; i++) {
            for (int j = 0; j < numberOfComponents; j++) {
                Assertions.assertEquals(matrixSequential[i][j], matrixParallel[i][j]);
            }
        }
    }

    /**
     * Tests calculation integrity across all 4 descriptor calculation methods for 1000 molecules.
     * Verifies that all input processing strategies (molecule parallelization, molecule batch parallelization,
     * SMILES strings batch parallelization and SMILES string parallelization) yield descriptor value matrices
     * identical to the sequential calculation with molecule parallelization. The correctness of the values themselves
     * is tested by the tests of the single descriptors.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    @Tag("SlowTest")
    void test_Integrity() throws Exception {
        int numberOfMolecules = 1000;
        IAtomContainer[] moleculesArray = createPropionicAcidMolecules(numberOfMolecules);
        String[] moleculeStringsArray = createPropionicAcidSmiles(numberOfMolecules);
        int startIndex = 0;
        Descriptor[] descriptors = Descriptor.getSpecifiedDescriptors(true, true, true, false);
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        boolean isParallelCalculation = false;
        int batchSize = 100;

        // reference: sequential calculation of the engine with molecule parallelization
        float[][] matrix2 = new float[numberOfMolecules][numberOfComponents];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
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
                CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
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
                CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
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
                Assertions.assertEquals(matrix2[i][j], matrix3[i][j]);
                Assertions.assertEquals(matrix2[i][j], matrix4[i][j]);
                Assertions.assertEquals(matrix2[i][j], matrix5[i][j]);
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

        // Test parallelization by single molecule parallelization
        float[][] matrixSequential = new float[numberOfMolecules][numberOfComponents];
        boolean isParallelCalculation = false;
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                        descriptors,
                        moleculesArray,
                        matrixSequential,
                        startIndex,
                        isParallelCalculation,
                        nanPositions
                )
        );

        float[][] matrixParallel = new float[numberOfMolecules][numberOfComponents];
        isParallelCalculation = true;
        List<int[]> nanPositionsParallel = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
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
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
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
     * Tests integrity for 3D coordinates descriptors across calculation methods: the batch parallelization yields
     * a matrix identical to the sequential calculation with molecule parallelization.
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

        // reference: sequential calculation of the engine with molecule parallelization
        float[][] matrix2 = new float[numberOfMolecules][numberOfComponents];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());
        Assertions.assertTrue(
                CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
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
                CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
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
                Assertions.assertEquals(matrix2[i][j], matrix3[i][j]);
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
        IAtomContainer methaneExplicit = DescriptorCalculator.createMoleculeWithExplicitHydrogens(methaneImplicit);

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
        IAtomContainer ethanolExplicit = DescriptorCalculator.createMoleculeWithExplicitHydrogens(ethanolImplicit);

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
        IAtomContainer benzeneClone = DescriptorCalculator.copyMolecule(benzene);

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

        IAtomContainer copy = DescriptorCalculator.copyMolecule(original);

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

        IAtomContainer copy = DescriptorCalculator.copyMolecule(original);
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

        IAtomContainer copy = DescriptorCalculator.copyMolecule(original);

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
     * Tests the configuration of the calculator: the default configuration, a non-default configuration and invalid
     * sizes.
     */
    @Test
    void test_Configuration() {
        Assertions.assertEquals(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, CALCULATOR.getFingerprintPoolSize());
        Assertions.assertEquals(4, CALCULATOR.getFingerprintPoolSize());
        Assertions.assertEquals(1024, CALCULATOR.getCircularFingerprintSize());
        Assertions.assertFalse(CALCULATOR.isCircularFingerprintCounts());

        DescriptorCalculator calculator = new DescriptorCalculator(6, 2048, true);
        Assertions.assertEquals(6, calculator.getFingerprintPoolSize());
        Assertions.assertEquals(2048, calculator.getCircularFingerprintSize());
        Assertions.assertTrue(calculator.isCircularFingerprintCounts());

        Assertions.assertThrows(IllegalArgumentException.class, () -> new DescriptorCalculator(0, 1024, false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> new DescriptorCalculator(-1, 1024, false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> new DescriptorCalculator(4, 0, false));
        Assertions.assertThrows(IllegalArgumentException.class, () -> new DescriptorCalculator(4, -1, false));
    }

    /**
     * Tests a calculator with a non-default circular fingerprint size (2048): the number of components of all circular
     * fingerprints is the configured size for this calculator, while the enum keeps the default size and other
     * descriptors are not affected; the calculated bits are equal to the bit fingerprints calculated directly with the
     * CDK class {@link CircularFingerprinter} with the same size, also for the descriptor after the circular fingerprint
     * (start indices); a matrix sized with the default size of the enum is too small.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CircularFingerprintSize() throws Exception {
        DescriptorCalculator calculator = new DescriptorCalculator(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, 2048, false);
        for (Descriptor descriptor : Descriptor.values()) {
            if (descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_")) {
                Assertions.assertEquals(2048, calculator.getNumberOfComponents(descriptor), descriptor.name());
                Assertions.assertEquals(1024, descriptor.getDescriptorComponentNumber(), descriptor.name());
                Assertions.assertEquals(1024, CALCULATOR.getNumberOfComponents(descriptor), descriptor.name());
            } else {
                Assertions.assertEquals(descriptor.getDescriptorComponentNumber(), calculator.getNumberOfComponents(descriptor), descriptor.name());
            }
        }

        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles("c1ccccc1CCc1ccccc1"); // 1,2-Diphenylethane CID: 7647
        DescriptorCalculator.setAromaticity(molecule, Aromaticity.Model.Daylight);
        IAtomContainer[] moleculesArray = new IAtomContainer[]{molecule};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, Descriptor.MOLECULAR_WEIGHT};
        int numberOfComponents = calculator.getNumberOfComponents(descriptors);
        Assertions.assertEquals(2049, numberOfComponents);
        Assertions.assertEquals(2049, calculator.getColumnNames(descriptors).length);
        Assertions.assertEquals("ECFP4_2047", calculator.getColumnNames(descriptors)[2047]);

        float[][] matrix = new float[1][numberOfComponents];
        List<int[]> nanPositions = new ArrayList<>();
        Assertions.assertTrue(calculator.setDescriptorsForMoleculesByMoleculeParallelization(descriptors, moleculesArray, matrix, 0, false, nanPositions));
        Assertions.assertTrue(nanPositions.isEmpty());
        IBitFingerprint fingerprint = new CircularFingerprinter(CircularFingerprinter.CLASS_ECFP4, 2048).getBitFingerprint(molecule);
        for (int i = 0; i < 2048; i++) {
            Assertions.assertEquals(fingerprint.get(i) ? 1.0f : 0.0f, matrix[0][i], "bit " + i);
        }
        Assertions.assertEquals(CALCULATOR.calculateDescriptor(Descriptor.MOLECULAR_WEIGHT, molecule)[0], matrix[0][2048]);
        Assertions.assertEquals(2048, calculator.calculateDescriptor(Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, molecule).length);

        float[][] tooSmallMatrix = new float[1][Descriptor.getNumberOfComponents(descriptors)];
        Assertions.assertThrows(IllegalArgumentException.class, () -> calculator.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, moleculesArray, tooSmallMatrix, 0, false, new ArrayList<>()));
    }

    /**
     * Tests getDescriptorAndComponentIndex() and getDescriptorAndComponentInfo() of the calculator: column indices
     * are mapped to the descriptor, the component index and the column name; with a non-default circular fingerprint
     * size, the columns after a circular fingerprint are shifted accordingly; invalid input.
     */
    @Test
    void test_getDescriptorAndComponentInfo() {
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.WIENER_NUMBER,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, Descriptor.ATOM_COUNT_C};
        Assertions.assertArrayEquals(new int[]{0, 0}, CALCULATOR.getDescriptorAndComponentIndex(descriptors, 0));
        Assertions.assertArrayEquals(new int[]{1, 1}, CALCULATOR.getDescriptorAndComponentIndex(descriptors, 2));
        Assertions.assertArrayEquals(new String[]{Descriptor.WIENER_NUMBER.getName(), "1", "WPOL"},
                CALCULATOR.getDescriptorAndComponentInfo(descriptors, 2));
        Assertions.assertArrayEquals(new String[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4.getName(), "5", "ECFP4_5"},
                CALCULATOR.getDescriptorAndComponentInfo(descriptors, 8));
        Assertions.assertArrayEquals(new String[]{Descriptor.ATOM_COUNT_C.getName(), "0", "nC"},
                CALCULATOR.getDescriptorAndComponentInfo(descriptors, 3 + 1024));
        // the default calculator gives the same indices as the enum
        for (int index : new int[]{0, 1, 2, 3, 500, 3 + 1023, 3 + 1024}) {
            Assertions.assertArrayEquals(Descriptor.getDescriptorAndComponentIndex(descriptors, index),
                    CALCULATOR.getDescriptorAndComponentIndex(descriptors, index), "index " + index);
        }

        DescriptorCalculator calculator = new DescriptorCalculator(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, 2048, false);
        Assertions.assertArrayEquals(new String[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4.getName(), "1024", "ECFP4_1024"},
                calculator.getDescriptorAndComponentInfo(descriptors, 3 + 1024));
        Assertions.assertArrayEquals(new String[]{Descriptor.ATOM_COUNT_C.getName(), "0", "nC"},
                calculator.getDescriptorAndComponentInfo(descriptors, 3 + 2048));

        Assertions.assertThrows(IllegalArgumentException.class, () -> CALCULATOR.getDescriptorAndComponentIndex(descriptors, -1));
        Assertions.assertThrows(IllegalArgumentException.class, () -> CALCULATOR.getDescriptorAndComponentIndex(descriptors, 3 + 1025));
        Assertions.assertThrows(IllegalArgumentException.class, () -> CALCULATOR.getDescriptorAndComponentIndex(new Descriptor[0], 0));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getDescriptorAndComponentIndex(null, 0));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getDescriptorAndComponentInfo(
                new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, null}, 0));
    }

    /**
     * Tests that the same molecule object must not occur more than once in the molecule array of both methods for
     * IAtomContainer arrays (the matrix is not changed), while copies of a molecule and identical SMILES strings are
     * allowed.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_DuplicateMolecules() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer molecule = smilesParser.parseSmiles("CC(=O)O"); //Acetic Acid CID: 176
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT};
        IAtomContainer[] duplicates = new IAtomContainer[]{molecule, molecule};
        float[][] matrix = new float[2][1];
        Assertions.assertThrows(IllegalArgumentException.class, () -> CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, duplicates, matrix, 0, true, new ArrayList<>()));
        Assertions.assertThrows(IllegalArgumentException.class, () -> CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                descriptors, duplicates, matrix, 0, 1, true, new ArrayList<>()));
        Assertions.assertArrayEquals(new float[2][1], matrix);

        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, new IAtomContainer[]{molecule, molecule.clone()}, matrix, 0, true, new ArrayList<>()));
        Assertions.assertEquals(matrix[0][0], matrix[1][0]);
        float[][] smilesMatrix = new float[2][1];
        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                descriptors, new String[]{"CC(=O)O", "CC(=O)O"}, smilesMatrix, 0, null, true, new ArrayList<>()));
        Assertions.assertArrayEquals(matrix, smilesMatrix);
    }

    /**
     * Comprehensive test verifying that all descriptor enum values return an empty float array (length 0)
     * when calculating descriptors for an empty molecule (0 atoms).
     */
    @Test
    void test_CalculateDescriptor_EmptyMolecule_AllDescriptors() {
        IAtomContainer emptyMolecule = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        Assertions.assertTrue(emptyMolecule.isEmpty());

        for (Descriptor descriptor : Descriptor.values()) {
            // Test calculateDescriptor(descriptor, emptyMolecule)
            float[] result1 = CALCULATOR.calculateDescriptor(descriptor, emptyMolecule);
            Assertions.assertNotNull(result1, "Result should not be null for " + descriptor);
            Assertions.assertEquals(0, result1.length, "Empty molecule should return 0-length array for " + descriptor);

            // Test calculateDescriptor(descriptor, emptyMolecule, electronDonationModel) with various models
            float[] result2 = CALCULATOR.calculateDescriptor(descriptor, emptyMolecule, ElectronDonation.daylight());
            Assertions.assertNotNull(result2, "Result should not be null for " + descriptor + " with Daylight model");
            Assertions.assertEquals(0, result2.length, "Empty molecule should return 0-length array for " + descriptor + " with Daylight model");

            float[] result3 = CALCULATOR.calculateDescriptor(descriptor, emptyMolecule, ElectronDonation.cdk());
            Assertions.assertNotNull(result3, "Result should not be null for " + descriptor + " with CDK model");
            Assertions.assertEquals(0, result3.length, "Empty molecule should return 0-length array for " + descriptor + " with CDK model");

            float[] result4 = CALCULATOR.calculateDescriptor(descriptor, emptyMolecule, ElectronDonation.piBonds());
            Assertions.assertNotNull(result4, "Result should not be null for " + descriptor + " with piBonds model");
            Assertions.assertEquals(0, result4.length, "Empty molecule should return 0-length array for " + descriptor + " with piBonds model");
        }
    }

    /**
     * Comprehensive test verifying that all descriptor enum values return an empty float array (length 0)
     * when calculating descriptors for empty and whitespace SMILES strings.
     */
    @Test
    void test_CalculateDescriptor_EmptyAndWhitespaceSmiles_AllDescriptors() {
        String[] emptySmilesVariants = new String[]{"", "   ", "\t", "\n", "\r\n", " \t \n \r "};

        for (Descriptor descriptor : Descriptor.values()) {
            for (String emptySmiles : emptySmilesVariants) {
                float[] result = CALCULATOR.calculateDescriptor(descriptor, emptySmiles, ElectronDonation.daylight());
                Assertions.assertNotNull(result, "Result should not be null for " + descriptor + " with SMILES: '" + emptySmiles + "'");
                Assertions.assertEquals(0, result.length, "Empty/whitespace SMILES should return 0-length array for " + descriptor + " with SMILES: '" + emptySmiles + "'");
            }
        }
    }

    /**
     * Tests that setAromaticity and createMoleculeWithExplicitHydrogens handle empty molecules gracefully.
     */
    @Test
    void test_SetAromaticity_And_ExplicitHydrogens_EmptyMolecule() throws Exception {
        IAtomContainer emptyMolecule = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        Assertions.assertTrue(emptyMolecule.isEmpty());

        // Should return early without throwing any exception
        DescriptorCalculator.setAromaticity(emptyMolecule, Aromaticity.Model.Daylight);
        Assertions.assertTrue(emptyMolecule.isEmpty());

        DescriptorCalculator.setAromaticity(emptyMolecule, Aromaticity.Model.CDK_1x);
        Assertions.assertTrue(emptyMolecule.isEmpty());

        // Should return empty copy/molecule without throwing any exception
        IAtomContainer explicitHydrogens = DescriptorCalculator.createMoleculeWithExplicitHydrogens(emptyMolecule);
        Assertions.assertNotNull(explicitHydrogens);
        Assertions.assertTrue(explicitHydrogens.isEmpty());
    }

    /**
     * Tests setDescriptorsForMoleculesByMoleculeParallelization with empty and mixed molecules in sequential and parallel modes.
     * Verifies that calculations complete safely without exceptions, valid molecules produce expected values,
     * and empty molecules are calculated consistently (e.g. 0.0 weight/count).
     */
    @Test
    void test_MoleculeParallelization_WithEmptyMolecules() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer ethanol = smilesParser.parseSmiles("CCO");
        IAtomContainer benzene = smilesParser.parseSmiles("c1ccccc1");
        IAtomContainer emptyMol1 = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);
        IAtomContainer emptyMol2 = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);

        IAtomContainer[] mixedMolecules = new IAtomContainer[]{ethanol, emptyMol1, benzene, emptyMol2};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT};
        int totalComponents = Descriptor.getNumberOfComponents(descriptors);

        for (boolean isParallel : new boolean[]{false, true}) {
            float[][] matrix = new float[mixedMolecules.length][totalComponents];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

            boolean success = CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                    descriptors, mixedMolecules, matrix, 0, isParallel, nanPositions
            );

            Assertions.assertTrue(success, "Calculation should succeed without NaNs (parallel=" + isParallel + ")");
            Assertions.assertEquals(0, nanPositions.size());

            // Valid molecules: ethanol (row 0) and benzene (row 2)
            Assertions.assertTrue(matrix[0][0] > 40.0f && matrix[0][0] < 50.0f, "Ethanol weight ~ 46");
            Assertions.assertEquals(9.0f, matrix[0][1], 0.01f, "Ethanol atom count with explicit H = 9");
            Assertions.assertTrue(matrix[2][0] > 70.0f && matrix[2][0] < 85.0f, "Benzene weight ~ 78");
            Assertions.assertEquals(12.0f, matrix[2][1], 0.01f, "Benzene atom count with explicit H = 12");

            // Empty molecules: emptyMol1 (row 1) and emptyMol2 (row 3) have 0.0 values
            Assertions.assertEquals(0.0f, matrix[1][0], 0.001f, "Empty molecule 1 weight = 0");
            Assertions.assertEquals(0.0f, matrix[1][1], 0.001f, "Empty molecule 1 atom count = 0");
            Assertions.assertEquals(0.0f, matrix[3][0], 0.001f, "Empty molecule 2 weight = 0");
            Assertions.assertEquals(0.0f, matrix[3][1], 0.001f, "Empty molecule 2 atom count = 0");
        }
    }

    /**
     * Tests setDescriptorsForMoleculesByBatchParallelization with empty and mixed molecules in sequential and parallel modes
     * across different batch sizes.
     */
    @Test
    void test_BatchParallelization_WithEmptyMolecules() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer ethanol = smilesParser.parseSmiles("CCO");
        IAtomContainer benzene = smilesParser.parseSmiles("c1ccccc1");
        IAtomContainer emptyMol = SilentChemObjectBuilder.getInstance().newInstance(IAtomContainer.class);

        IAtomContainer[] molecules = new IAtomContainer[]{ethanol, emptyMol, benzene};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT};
        int totalComponents = Descriptor.getNumberOfComponents(descriptors);

        for (boolean isParallel : new boolean[]{false, true}) {
            for (int batchSize : new int[]{1, 2, 3}) {
                float[][] matrix = new float[molecules.length][totalComponents];
                List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

                boolean success = CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                        descriptors, molecules, matrix, 0, batchSize, isParallel, nanPositions
                );

                Assertions.assertTrue(success, "Calculation should succeed (batchSize=" + batchSize + ", parallel=" + isParallel + ")");
                Assertions.assertTrue(matrix[0][0] > 40.0f, "Ethanol weight should be positive");
                Assertions.assertEquals(0.0f, matrix[1][0], 0.001f, "Empty molecule weight should be 0.0");
                Assertions.assertEquals(0.0f, matrix[1][1], 0.001f, "Empty molecule atom count should be 0.0");
                Assertions.assertTrue(matrix[2][0] > 70.0f, "Benzene weight should be positive");
            }
        }
    }

    /**
     * Tests setDescriptorsForMoleculesBySmilesStringParallelization with empty and whitespace SMILES in sequential and parallel modes.
     */
    @Test
    void test_SmilesParallelization_WithEmptyAndWhitespaceSmiles() throws Exception {
        String[] smilesArray = new String[]{"CCO", "", "c1ccccc1", "   ", "\t\n", "CC(=O)O"};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT};
        int totalComponents = Descriptor.getNumberOfComponents(descriptors);

        for (boolean isParallel : new boolean[]{false, true}) {
            float[][] matrix = new float[smilesArray.length][totalComponents];
            List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

            boolean success = CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                    descriptors, smilesArray, matrix, 0, ElectronDonation.daylight(), isParallel, nanPositions
            );

            Assertions.assertTrue(success, "Calculation with empty/whitespace SMILES should succeed without NaNs (parallel=" + isParallel + ")");

            // Valid SMILES rows: CCO (row 0), Benzene (row 2), Acetic acid (row 5)
            Assertions.assertTrue(matrix[0][0] > 40.0f, "CCO weight > 40");
            Assertions.assertEquals(9.0f, matrix[0][1], 0.01f, "CCO atom count with explicit H = 9");
            Assertions.assertTrue(matrix[2][0] > 70.0f, "Benzene weight > 70");
            Assertions.assertEquals(12.0f, matrix[2][1], 0.01f, "Benzene atom count with explicit H = 12");
            Assertions.assertTrue(matrix[5][0] > 55.0f, "Acetic acid weight > 55");
            Assertions.assertEquals(8.0f, matrix[5][1], 0.01f, "Acetic acid atom count with explicit H = 8");

            // Empty/whitespace SMILES rows: empty (row 1), spaces (row 3), whitespace (row 4)
            Assertions.assertEquals(0.0f, matrix[1][0], 0.001f, "Empty SMILES weight = 0");
            Assertions.assertEquals(0.0f, matrix[1][1], 0.001f, "Empty SMILES atom count = 0");
            Assertions.assertEquals(0.0f, matrix[3][0], 0.001f, "Spaces SMILES weight = 0");
            Assertions.assertEquals(0.0f, matrix[3][1], 0.001f, "Spaces SMILES atom count = 0");
            Assertions.assertEquals(0.0f, matrix[4][0], 0.001f, "Whitespace SMILES weight = 0");
            Assertions.assertEquals(0.0f, matrix[4][1], 0.001f, "Whitespace SMILES atom count = 0");
        }
    }

    /**
     * Tests setDescriptorsForMoleculeBySmilesStringsBatchParallelization with empty and whitespace SMILES.
     */
    @Test
    void test_SmilesBatchParallelization_WithEmptyAndWhitespaceSmiles() throws Exception {
        String[] smilesArray = new String[]{"CCO", "", "c1ccccc1", "   "};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.ATOM_COUNT};
        int totalComponents = Descriptor.getNumberOfComponents(descriptors);

        for (boolean isParallel : new boolean[]{false, true}) {
            for (int batchSize : new int[]{1, 2, 4}) {
                float[][] matrix = new float[smilesArray.length][totalComponents];
                List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

                boolean success = CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                        descriptors, smilesArray, matrix, 0, batchSize, ElectronDonation.daylight(), isParallel, nanPositions
                );

                Assertions.assertTrue(success, "Calculation should succeed (batchSize=" + batchSize + ", parallel=" + isParallel + ")");
                Assertions.assertTrue(matrix[0][0] > 40.0f, "CCO weight should be positive");
                Assertions.assertEquals(0.0f, matrix[1][0], 0.001f, "Empty SMILES should be 0.0");
                Assertions.assertEquals(0.0f, matrix[1][1], 0.001f, "Empty SMILES atom count should be 0.0");
                Assertions.assertTrue(matrix[2][0] > 70.0f, "Benzene weight should be positive");
                Assertions.assertEquals(0.0f, matrix[3][0], 0.001f, "Spaces SMILES should be 0.0");
                Assertions.assertEquals(0.0f, matrix[3][1], 0.001f, "Spaces SMILES atom count should be 0.0");
            }
        }
    }

    /**
     * Tests that empty descriptor arrays, empty molecule arrays, and empty SMILES arrays
     * abort calculation gracefully and return false.
     */
    @Test
    void test_EmptyArraysValidation() throws Exception {
        Descriptor[] validDescriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT};
        Descriptor[] emptyDescriptors = new Descriptor[0];

        IAtomContainer[] validMolecules = new IAtomContainer[]{
                (new SmilesParser(SilentChemObjectBuilder.getInstance())).parseSmiles("CCO")
        };
        IAtomContainer[] emptyMolecules = new IAtomContainer[0];

        String[] validSmiles = new String[]{"CCO"};
        String[] emptySmiles = new String[0];

        float[][] matrix = new float[1][1];
        List<int[]> nanPositions = Collections.synchronizedList(new LinkedList<>());

        // Empty descriptors array
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                emptyDescriptors, validMolecules, matrix, 0, false, nanPositions
        ));
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                emptyDescriptors, validMolecules, matrix, 0, 1, false, nanPositions
        ));
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                emptyDescriptors, validSmiles, matrix, 0, ElectronDonation.daylight(), false, nanPositions
        ));
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                emptyDescriptors, validSmiles, matrix, 0, 1, ElectronDonation.daylight(), false, nanPositions
        ));

        // Empty molecule array
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                validDescriptors, emptyMolecules, matrix, 0, false, nanPositions
        ));
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                validDescriptors, emptyMolecules, matrix, 0, 1, false, nanPositions
        ));

        // Empty SMILES array
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                validDescriptors, emptySmiles, matrix, 0, ElectronDonation.daylight(), false, nanPositions
        ));
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculeBySmilesStringsBatchParallelization(
                validDescriptors, emptySmiles, matrix, 0, 1, ElectronDonation.daylight(), false, nanPositions
        ));
    }

    /**
     * Tests that NaN positions are tracked in any (non-synchronized) list, also in parallel mode, and are added
     * sorted by molecule index and component index, including the start index offset.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_NaNPositions_PlainListSortedWithStartIndex() throws Exception {
        // HYBRIDIZATION_RATIO is NaN for "[Na+].[Cl-]" and "[Xe]"
        String[] pattern = new String[]{"CCO", "[Na+].[Cl-]", "CC", "[Xe]"};
        String[] smiles = new String[200];
        for (int i = 0; i < smiles.length; i++) {
            smiles[i] = pattern[i % pattern.length];
        }
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.HYBRIDIZATION_RATIO};
        int startIndex = 1;
        int nanColumn = startIndex + Descriptor.MOLECULAR_WEIGHT.getDescriptorComponentNumber();
        List<int[]> expected = new ArrayList<>();
        for (int i = 0; i < smiles.length; i++) {
            if (i % 2 == 1) {
                expected.add(new int[]{i, nanColumn});
            }
        }
        for (boolean isParallel : new boolean[]{false, true}) {
            float[][] matrix = new float[smiles.length][startIndex + Descriptor.getNumberOfComponents(descriptors)];
            List<int[]> nanPositions = new ArrayList<>();
            boolean success = CALCULATOR.setDescriptorsForMoleculesBySmilesStringParallelization(
                    descriptors, smiles, matrix, startIndex, null, isParallel, nanPositions);
            Assertions.assertFalse(success, "NaN values expected (parallel=" + isParallel + ")");
            Assertions.assertEquals(expected.size(), nanPositions.size(), "parallel=" + isParallel);
            for (int i = 0; i < expected.size(); i++) {
                Assertions.assertArrayEquals(expected.get(i), nanPositions.get(i), "position " + i + ", parallel=" + isParallel);
                Assertions.assertTrue(Float.isNaN(matrix[expected.get(i)[0]][expected.get(i)[1]]));
            }
        }
    }

    /**
     * Tests that the calculation recipe of every descriptor matches its type, i.e. fingerprints are calculated with a
     * pool of fingerprinters and all other descriptors with a CDK molecular descriptor.
     */
    @Test
    void test_CalculationFactoryMatchesDescriptorType() {
        for (Descriptor descriptor : Descriptor.values()) {
            DescriptorCalculation calculation = descriptor.getCalculationFactory().create(
                    new CalculationSettings(1, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, false));
            Assertions.assertEquals(descriptor.isFingerprint(), calculation instanceof FingerprintCalculation, descriptor.name());
            Assertions.assertEquals(!descriptor.isFingerprint(), calculation instanceof MolecularCalculation, descriptor.name());
            // in count mode, only circular fingerprints are calculated as counts
            DescriptorCalculation countCalculation = descriptor.getCalculationFactory().create(
                    new CalculationSettings(1, Descriptor.CIRCULAR_FINGERPRINT_DEFAULT_SIZE, true));
            Assertions.assertEquals(descriptor.name().startsWith("CIRCULAR_FINGERPRINTER_"),
                    countCalculation instanceof FingerprintCalculation && ((FingerprintCalculation) countCalculation).isCounts(),
                    descriptor.name());
        }
    }

    /**
     * Tests that the 3D precondition is checked for all molecules, not only for the first one: if any molecule lacks
     * 3D coordinates, the calculation is aborted (false is returned and the matrix is not changed).
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_3DPrecondition_CheckedForAllMolecules() throws Exception {
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer with3D = smilesParser.parseSmiles("CC");
        with3D.getAtom(0).setPoint3d(new Point3d(0.0, 0.0, 0.0));
        with3D.getAtom(1).setPoint3d(new Point3d(1.54, 0.0, 0.0));
        IAtomContainer without3D = smilesParser.parseSmiles("CC");
        IAtomContainer[] molecules = new IAtomContainer[]{with3D, with3D.clone(), without3D};
        Descriptor[] descriptors = new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.MOMENT_OF_INERTIA};

        float[][] matrix = new float[molecules.length][Descriptor.getNumberOfComponents(descriptors)];
        List<int[]> nanPositions = new ArrayList<>();
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, molecules, matrix, 0, true, nanPositions));
        float[][] batchMatrix = new float[molecules.length][Descriptor.getNumberOfComponents(descriptors)];
        Assertions.assertFalse(CALCULATOR.setDescriptorsForMoleculesByBatchParallelization(
                descriptors, molecules, batchMatrix, 0, 2, true, nanPositions));

        Assertions.assertTrue(nanPositions.isEmpty());
        for (int i = 0; i < molecules.length; i++) {
            Assertions.assertArrayEquals(new float[matrix[i].length], matrix[i], "calculation must not have started");
            Assertions.assertArrayEquals(new float[batchMatrix[i].length], batchMatrix[i], "calculation must not have started");
        }
    }

    /**
     * Tests the folded count mode of the circular fingerprints (see
     * {@link DescriptorCalculator#DescriptorCalculator(int, int, boolean)}):
     * <ul>
     *     <li>the folding of the hash codes is checked against reference values of the CDK test class
     *         CircularFingerprinterTest: folding the reference hash codes and counts of the count fingerprint of
     *         propanamide ("CCC(=O)N", ECFP6) must yield exactly the reference bits of its bit fingerprint
     *         (1024 bits), and the engine must calculate exactly these counts;</li>
     *     <li>for all ECFP and FCFP fingerprints, a count is greater than 0 if and only if the corresponding bit of the
     *         bit fingerprint is set;</li>
     *     <li>all counts are equal to the count fingerprints calculated directly with the CDK class
     *         {@link CircularFingerprinter}, i.e. without the engine, and folded;</li>
     *     <li>counts greater than 1 occur (1,2-diphenylethane), and the number of components does not change.</li>
     * </ul>
     * Note: This test checks that the engine calculates the fingerprints correctly, NOT that the fingerprints
     * themselves are scientifically correct; the latter is the responsibility of the CDK fingerprint tests.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    void test_CircularFingerprintCounts() throws Exception {
        DescriptorCalculator countCalculator = new DescriptorCalculator(DescriptorCalculator.DEFAULT_FINGERPRINT_POOL_SIZE, 1024, true);
        int size = countCalculator.getCircularFingerprintSize();
        Assertions.assertEquals(1024, size, "reference values require the default size of 1024 bits");
        Assertions.assertTrue(countCalculator.isCircularFingerprintCounts());
        Assertions.assertFalse(CALCULATOR.isCircularFingerprintCounts(), "bit fingerprints are the default");

        // reference values from CircularFingerprinterTest (propanamide, ECFP6, 1024 bits)
        int[] referenceHashCodesAndCounts = new int[]{-414937772, 1, -1027418143, 1, 1627608083, 1, -868007456, 1,
                -1006701866, 1, -1059145289, 1, -801752141, 1, 790592664, 1, -289109509, 1, -1650154758, 1, 1286833445, 1};
        int[] referenceBits = new int[]{19, 152, 293, 340, 439, 480, 507, 726, 762, 947, 993};
        float[] expectedPropanamideEcfp6Counts = new float[size];
        BitSet foldedBits = new BitSet();
        for (int i = 0; i < referenceHashCodesAndCounts.length; i += 2) {
            int position = FingerprintCalculation.fold(referenceHashCodesAndCounts[i], size);
            expectedPropanamideEcfp6Counts[position] += referenceHashCodesAndCounts[i + 1];
            foldedBits.set(position);
        }
        BitSet expectedBits = new BitSet();
        for (int bit : referenceBits) {
            expectedBits.set(bit);
        }
        Assertions.assertEquals(expectedBits, foldedBits, "folding of the reference count fingerprint must yield the reference bits");

        Descriptor[] descriptors = new Descriptor[]{
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_0, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4, Descriptor.CIRCULAR_FINGERPRINTER_ECFP_6,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_0, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_2,
                Descriptor.CIRCULAR_FINGERPRINTER_FCFP_4, Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6
        };
        int[] cdkClasses = new int[]{
                CircularFingerprinter.CLASS_ECFP0, CircularFingerprinter.CLASS_ECFP2,
                CircularFingerprinter.CLASS_ECFP4, CircularFingerprinter.CLASS_ECFP6,
                CircularFingerprinter.CLASS_FCFP0, CircularFingerprinter.CLASS_FCFP2,
                CircularFingerprinter.CLASS_FCFP4, CircularFingerprinter.CLASS_FCFP6
        };
        SmilesParser smilesParser = new SmilesParser(SilentChemObjectBuilder.getInstance());
        IAtomContainer diphenylethane = smilesParser.parseSmiles("c1ccccc1CCc1ccccc1"); // 1,2-Diphenylethane CID: 7647
        DescriptorCalculator.setAromaticity(diphenylethane, Aromaticity.Model.Daylight);
        IAtomContainer propanamide = smilesParser.parseSmiles("CCC(=O)N"); // Propanamide CID: 6578, as in the CDK test
        IAtomContainer[] moleculesArray = new IAtomContainer[]{diphenylethane, propanamide};
        int numberOfComponents = Descriptor.getNumberOfComponents(descriptors);
        Assertions.assertEquals(descriptors.length * size, numberOfComponents);

        // reference: the count fingerprints calculated directly with the CDK class, without the engine, and folded
        float[][] expectedCounts = new float[moleculesArray.length][numberOfComponents];
        for (int m = 0; m < moleculesArray.length; m++) {
            for (int d = 0; d < descriptors.length; d++) {
                ICountFingerprint countFingerprint = new CircularFingerprinter(cdkClasses[d], size).getCountFingerprint(moleculesArray[m]);
                for (int i = 0; i < countFingerprint.numOfPopulatedbins(); i++) {
                    expectedCounts[m][d * size + FingerprintCalculation.fold(countFingerprint.getHash(i), size)] += countFingerprint.getCount(i);
                }
            }
        }

        float[][] bitMatrix = new float[moleculesArray.length][numberOfComponents];
        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, moleculesArray, bitMatrix, 0, false, new ArrayList<>()));
        Assertions.assertEquals(numberOfComponents, countCalculator.getNumberOfComponents(descriptors), "number of components must not change");
        for (boolean isParallelCalculation : new boolean[]{false, true}) {
            float[][] countMatrix = new float[moleculesArray.length][numberOfComponents];
            List<int[]> nanPositions = new ArrayList<>();
            Assertions.assertTrue(countCalculator.setDescriptorsForMoleculesByMoleculeParallelization(
                    descriptors, moleculesArray, countMatrix, 0, isParallelCalculation, nanPositions));
            Assertions.assertTrue(nanPositions.isEmpty());
            boolean hasCountGreaterThanOne = false;
            for (int m = 0; m < moleculesArray.length; m++) {
                Assertions.assertArrayEquals(expectedCounts[m], countMatrix[m], "molecule " + m + ", parallel=" + isParallelCalculation);
                for (int i = 0; i < numberOfComponents; i++) {
                    Assertions.assertEquals(bitMatrix[m][i] == 1.0f, countMatrix[m][i] > 0.0f, "molecule " + m + ", component " + i);
                    hasCountGreaterThanOne |= countMatrix[m][i] > 1.0f;
                }
            }
            Assertions.assertTrue(hasCountGreaterThanOne, "counts greater than 1 expected for 1,2-diphenylethane");
            // reference counts of propanamide (row 1), ECFP6 is the fourth fingerprint
            Assertions.assertArrayEquals(expectedPropanamideEcfp6Counts, Arrays.copyOfRange(countMatrix[1], 3 * size, 4 * size));
        }
        // the default calculator is not affected by the count calculator
        Assertions.assertFalse(CALCULATOR.isCircularFingerprintCounts());
        float[][] bitMatrixAgain = new float[moleculesArray.length][numberOfComponents];
        Assertions.assertTrue(CALCULATOR.setDescriptorsForMoleculesByMoleculeParallelization(
                descriptors, moleculesArray, bitMatrixAgain, 0, false, new ArrayList<>()));
        Assertions.assertArrayEquals(bitMatrix, bitMatrixAgain);
    }

    /**
     * Tests getColumnNames(): for every descriptor, the number of names is equal to its number of components and no
     * name is empty; the names of all descriptors are unique (the amino acid counts get the prefix "AA." because their
     * CDK names like "nC" collide with the atom counts); for several descriptors, the names are returned in the order
     * of the matrix columns (example with the CDK descriptor names "MW", "WPATH", "WPOL" and "nC"); fingerprints get
     * generated names (prefix and bit index); empty and null input.
     */
    @Test
    void test_getColumnNames() {
        for (Descriptor descriptor : Descriptor.values()) {
            Descriptor[] singleDescriptor = new Descriptor[]{descriptor};
            String[] names = CALCULATOR.getColumnNames(singleDescriptor);
            Assertions.assertEquals(Descriptor.getNumberOfComponents(singleDescriptor), names.length, descriptor.name());
            for (String name : names) {
                Assertions.assertNotNull(name, descriptor.name());
                Assertions.assertFalse(name.trim().isEmpty(), descriptor.name());
            }
        }
        Descriptor[] allDescriptors = Descriptor.values();
        String[] allNames = CALCULATOR.getColumnNames(allDescriptors);
        Assertions.assertEquals(Descriptor.getNumberOfComponents(allDescriptors), allNames.length);
        Set<String> uniqueNames = new HashSet<>();
        for (String name : allNames) {
            Assertions.assertTrue(uniqueNames.add(name), "Duplicate column name: " + name);
        }

        // amino acid counts: prefixed CDK names, atom counts unchanged
        String[] aminoAcidNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.AMINO_ACID_COUNT});
        Assertions.assertEquals("AA.nA", aminoAcidNames[0]);
        Assertions.assertTrue(Arrays.asList(aminoAcidNames).contains("AA.nC"));
        Assertions.assertArrayEquals(new String[]{"nC"}, CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.ATOM_COUNT_C}));

        // names in the order of the matrix columns
        Assertions.assertArrayEquals(new String[]{"MW", "WPATH", "WPOL", "nC"}, CALCULATOR.getColumnNames(
                new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, Descriptor.WIENER_NUMBER, Descriptor.ATOM_COUNT_C}));
        Assertions.assertArrayEquals(new String[]{"nC", "MW"}, CALCULATOR.getColumnNames(
                new Descriptor[]{Descriptor.ATOM_COUNT_C, Descriptor.MOLECULAR_WEIGHT}));

        // fingerprints: generated names
        String[] pubChemNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.PUBCHEM_FINGERPRINTER});
        Assertions.assertEquals("PubChem_0", pubChemNames[0]);
        Assertions.assertEquals("PubChem_880", pubChemNames[880]);
        String[] maccsNames = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.MACCS_FINGERPRINTER});
        Assertions.assertEquals("MACCS_0", maccsNames[0]);
        Assertions.assertEquals("MACCS_165", maccsNames[165]);
        int circularSize = CALCULATOR.getCircularFingerprintSize();
        String[] ecfp4Names = CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_ECFP_4});
        Assertions.assertEquals(circularSize, ecfp4Names.length);
        Assertions.assertEquals("ECFP4_0", ecfp4Names[0]);
        Assertions.assertEquals("ECFP4_" + (circularSize - 1), ecfp4Names[circularSize - 1]);
        Assertions.assertEquals("FCFP6_0", CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.CIRCULAR_FINGERPRINTER_FCFP_6})[0]);

        // empty and null input
        Assertions.assertEquals(0, CALCULATOR.getColumnNames(new Descriptor[0]).length);
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getColumnNames(null));
        Assertions.assertThrows(NullPointerException.class, () -> CALCULATOR.getColumnNames(new Descriptor[]{Descriptor.MOLECULAR_WEIGHT, null}));
    }
}
