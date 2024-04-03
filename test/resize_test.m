classdef resize_test < matlab.unittest.TestCase

    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
    end

    methods (Test)
        function test_pad_zeros(testCase)
            % pad_to_size example 1
            a = [1, 2, 3; 4, 5, 6];

            actual = pad_to_size(a, [5, 5]);
            expected = [
                0     0     0     0     0
                0     0     0     0     0
                0     1     2     3     0
                0     4     5     6     0
                0     0     0     0     0
                ];

            testCase.assertEqual(actual, expected);
        end

        function test_pad_val(testCase)
            % pad_to_size example 2
            a = [1, 2, 3; 4, 5, 6];

            actual = pad_to_size(a, [4, 6], -10);
            expected = [
                -10   -10   -10   -10   -10   -10
                -10   -10     1     2     3   -10
                -10   -10     4     5     6   -10
                -10   -10   -10   -10   -10   -10
                ];

            testCase.assertEqual(actual, expected);
        end

        function test_crop(testCase)
            % crop_to_size example
            cropIm = crop_to_size(zeros([8, 8, 8, 8]), [3, 4, 5]);

            testCase.assertEqual(size(cropIm), [3, 4, 5, 8]);
        end

        function test_inverse(testCase)
            % Test that crop_to_size exactly undoes the changes made by
            % pad_to_size
            for iter = 1:100
                sz = randi(5, 1, 5);
                szPad = randi(5, 1, 5) + 4;  % smallest of this line is equal to largest of previous

                initial = randn(sz);
                padded = pad_to_size(initial, szPad);
                cropped = crop_to_size(padded, sz);

                testCase.assertEqual(cropped, initial);
            end
        end
    end

end