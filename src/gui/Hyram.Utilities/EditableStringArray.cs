// Copyright 2016 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
// Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
// 
// This file is part of HyRAM (Hydrogen Risk Assessment Models).
// 
// HyRAM is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// HyRAM is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with HyRAM.  If not, see <https://www.gnu.org/licenses/>.

namespace SandiaNationalLaboratories.Hyram
{
    public enum ArrayStringConversionOption
    {
        AppendCarriageReturn,
        NoModifications,
        AppendCrlf
    }

    /// <summary>
    ///     Encapsulates a string array that can be resized, searched and manipulated without
    ///     having to call external functions or assign the result of a method call on one array
    ///     to a new one.  The Data property contains the array.
    /// </summary>
    public class EditableStringArray
    {

        /// <summary>
        ///     Create a new instance by referencing a preexisting string array.
        /// </summary>
        /// <param name="startingData">Data to reference.</param>
        public EditableStringArray(string[] startingData)
        {
            Data = startingData;
        }

        /// <summary>
        ///     Default constructor.  A new array with zero elements is created.
        /// </summary>
        public EditableStringArray()
        {
        }

        public string this[int index]
        {
            get => Data[index];
            set => Data[index] = value;
        }

        /// <summary>
        ///     The string array that this class operates upon.
        /// </summary>
        public string[] Data { get; set; } = new string[0];

        /// Insert an array into the Data Array.
        /// </summary>
        /// <param name="linesToInsert">An array containing the data to insert.</param>
        /// <param name="beforeIndex">Point at which to make insertion.</param>
        public void Insert(string[] linesToInsert, int beforeIndex)
        {
            Data = StringFunctions.InsertDataIntoArray(Data, linesToInsert, beforeIndex);
        }

        public string CombineToString(ArrayStringConversionOption conversionOption)
        {
            var result = "";
            var appendValue = "";

            if (conversionOption == ArrayStringConversionOption.AppendCarriageReturn)
                appendValue = "\r";
            else if (conversionOption == ArrayStringConversionOption.AppendCrlf) appendValue = "\r\n";

            foreach (var thisLine in Data)
                if (result != "")
                    result += appendValue + thisLine;
                else
                    result += thisLine;

            return result;
        }


        /// <summary>
        ///     Inserts one line into Data Array.
        /// </summary>
        /// <param name="lineToInsert">Line of data to insert.</param>
        /// <param name="beforeIndex">Index to insert lines before.</param>
        public void Insert(string lineToInsert, int beforeIndex)
        {
            var insertLines = new string[1];
            insertLines[0] = lineToInsert;
            Insert(insertLines, beforeIndex);
        }

        /// <summary>
        ///     Add a new line to the end of the Data Array.
        /// </summary>
        /// <param name="elementToAppend">New line to append.</param>
        public void Append(string elementToAppend)
        {
            var elem = new string[1];
            elem[0] = elementToAppend;
            Append(elem);
        }

        /// <summary>
        ///     Add an array of lines to the end of the Data Array.
        /// </summary>
        /// <param name="dataToAppend">The array of lines to append.</param>
        public void Append(string[] dataToAppend)
        {
            var newData = new string[dataToAppend.Length + Data.Length];
            Data.CopyTo(newData, 0);
            dataToAppend.CopyTo(newData, Data.Length);
            Data = newData;
        }
    }
}