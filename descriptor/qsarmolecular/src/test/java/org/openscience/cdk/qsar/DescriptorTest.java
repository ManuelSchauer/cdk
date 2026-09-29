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

/**
 * Test class for the {@link Descriptor} enum, i.e. the descriptor metadata and the static query methods of
 * the descriptor catalog. The calculation of descriptor values is tested in {@link DescriptorCalculatorTest}.
 *
 * @author Achim Zielesny
 * @author Jonas Schaub
 * @author Manuel Schauer
 */
class DescriptorTest {
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
