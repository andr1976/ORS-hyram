﻿// Copyright 2016 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
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

using System;
using System.Collections.Generic;
using System.Diagnostics;
using Python.Runtime;

namespace SandiaNationalLaboratories.Hyram
{
    /// <summary>
    ///     Results for single leak scenario of QRA analysis
    /// </summary>
    [Serializable]
    public class LeakResult
    {
        public string LeakSize;
        public double CompressorLeakFreq;
        public double CylinderLeakFreq;
        public double ExplosAvgEvents;
        public double ExplosionPllContrib;
        public double ExtraComp1LeakFreq;
        public double ExtraComp2LeakFreq;
        public double FilterLeakFreq;
        public double FlangeLeakFreq;
        public double HoseLeakFreq;
        public double InstrumentLeakFreq;
        public double JointLeakFreq;
        public double PipeLeakFreq;
        public double ValveLeakFreq;

        // -1 if not used
        public double H2ReleaseOverride;

        // Accident and shutdown failures for 100% leak size
        public double ProbOverpressureRupture;
        public double FreqOverpressureRupture;
        public double ProbDriveoffs;
        public double FreqDriveoffs;
        public double ProbNozzleRelease;
        public double FreqNozzleRelease;
        public double ProbMValveFtc;
        public double FreqMValveFtc;
        public double ProbSolValvesFtc;
        public double FreqSolValvesFtc;
        public double TotalFreqOtherFailures;

        public double ProbExplosion;
        public double ProbJetfire;
        public double ProbNoIgnition;
        public double ProbShutdown;

        public double JetfireAvgEvents;
        public double JetfirePllContrib;
        public double NoIgnAvgEvents;
        public double ShutdownAvgEvents;

        // Flag to check if override given. If ManualOverride is true, then only TotalProb will be shown.
        public double VehicleFailureOverride;

        public string GetLeakSizeString()
        {
            return LeakSize + "% Release";
        }

        public override string ToString()
        {
            var summary = $@"Leak result {LeakSize}:
                P jetfire: {ProbJetfire}
                P explos: {ProbExplosion}
                P no ign: {ProbNoIgnition}
                P shutd: {ProbShutdown}
                ";
            return summary;
        }
    }

    /// <summary>
    ///     Container for QRA analysis results returned by Python call
    /// </summary>
    [Serializable]
    public class QraResult
    {
        public double Air;
        public double Far;
        public List<LeakResult> LeakResults;
        public string[] PositionPlotFilenames;
        public double TotalPll;

        public QraResult(PyObject pyResult)
        {
            // Convert PyObject attrs to double with double-cast
            Air = (double) (dynamic) pyResult["air"];
            TotalPll = (double) (dynamic) pyResult["total_pll"];
            Far = (double) (dynamic) pyResult["far"];
            PositionPlotFilenames = (string[]) (dynamic) pyResult["plot_files"];

            // Parse scenario data for each leak size into objects
            dynamic leakResultData = pyResult["leak_results"];
            LeakResults = new List<LeakResult>();

            for (var i = 0; i < 5; i++)
            {
                var res = leakResultData[i];
                var leakSize = ((double) res["leak_size"]).ToString("000.00");
                var probShutdown = (double) res["p_shutdown"];
                var probJetfire = (double) res["p_jetfire"];
                var probExplosion = (double) res["p_explos"];
                var probNoIgnition = (double) res["p_no_ign"];

                var jetfireAvgEvents = (double) res["jetfire_avg_events"];
                var explosAvgEvents = (double) res["explos_avg_events"];
                var shutdownAvgEvents = (double) res["shutdown_avg_events"];
                var noIgnAvgEvents = (double) res["no_ign_avg_events"];

                var explosionPllContrib = (double) res["explos_pll_contrib"];
                var jetfirePllContrib = (double) res["jetfire_pll_contrib"];
                var h2ReleaseOverride = (double) res["release_freq_override"];

                var compressorLeakFreq = (double) res["component_leak_freqs"]["compressor"];
                var cylinderLeakFreq = (double) res["component_leak_freqs"]["cylinder"];
                var valveLeakFreq = (double) res["component_leak_freqs"]["valve"];
                var instrumentLeakFreq = (double) res["component_leak_freqs"]["instrument"];
                var jointLeakFreq = (double) res["component_leak_freqs"]["joint"];
                var hoseLeakFreq = (double) res["component_leak_freqs"]["hose"];
                var pipeLeakFreq = (double) res["component_leak_freqs"]["pipe"];
                var filterLeakFreq = (double) res["component_leak_freqs"]["filter"];
                var flangeLeakFreq = (double) res["component_leak_freqs"]["flange"];
                var extraComp1LeakFreq = (double) res["component_leak_freqs"]["extra1"];
                var extraComp2LeakFreq = (double) res["component_leak_freqs"]["extra2"];

                var nextLeakRes = new LeakResult
                {
                    LeakSize = leakSize,
                    ProbJetfire = probJetfire,
                    JetfireAvgEvents = jetfireAvgEvents,
                    JetfirePllContrib = jetfirePllContrib,
                    ProbExplosion = probExplosion,
                    ExplosAvgEvents = explosAvgEvents,
                    ShutdownAvgEvents = shutdownAvgEvents,
                    NoIgnAvgEvents = noIgnAvgEvents,
                    ExplosionPllContrib = explosionPllContrib,
                    ProbNoIgnition = probNoIgnition,
                    ProbShutdown = probShutdown,
                    H2ReleaseOverride = h2ReleaseOverride,

                    CompressorLeakFreq = compressorLeakFreq,
                    CylinderLeakFreq = cylinderLeakFreq,
                    ValveLeakFreq = valveLeakFreq,
                    InstrumentLeakFreq = instrumentLeakFreq,
                    JointLeakFreq = jointLeakFreq,
                    HoseLeakFreq = hoseLeakFreq,
                    PipeLeakFreq = pipeLeakFreq,
                    FilterLeakFreq = filterLeakFreq,
                    FlangeLeakFreq = flangeLeakFreq,
                    ExtraComp1LeakFreq = extraComp1LeakFreq,
                    ExtraComp2LeakFreq = extraComp2LeakFreq
                };

                if (i == 4)
                {
                    // Grab shutdown/accident failure data for 100% leak size
                    nextLeakRes.VehicleFailureOverride = (double) res["fueling_fail_freq_override"];
                    if (nextLeakRes.VehicleFailureOverride == -1.0)
                    {
                        nextLeakRes.ProbOverpressureRupture = (double) res["p_overp_rupture"];
                        nextLeakRes.FreqOverpressureRupture = (double) res["f_overp_rupture"];

                        nextLeakRes.ProbDriveoffs = (double) res["p_driveoff"];
                        nextLeakRes.FreqDriveoffs = (double) res["f_driveoff"];

                        nextLeakRes.ProbSolValvesFtc = (double) res["p_sol_valves_ftc"];
                        nextLeakRes.FreqSolValvesFtc = (double) res["f_sol_valves_ftc"];

                        nextLeakRes.ProbMValveFtc = (double) res["p_mvalve_ftc"];
                        nextLeakRes.FreqMValveFtc = (double) res["f_mvalve_ftc"];

                        nextLeakRes.ProbNozzleRelease = (double) res["p_nozzle_release"];
                        nextLeakRes.FreqNozzleRelease = (double) res["f_nozzle_release"];

                        nextLeakRes.TotalFreqOtherFailures = (double) res["fueling_fail_freq"];
                    }
                    else
                    {
                        nextLeakRes.TotalFreqOtherFailures = nextLeakRes.VehicleFailureOverride;
                    }
                }

                Trace.TraceInformation(nextLeakRes.ToString());
                LeakResults.Add(nextLeakRes);
            }
        }
    }
}