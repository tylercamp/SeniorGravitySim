using System;
using UnityEngine;
using System.Collections;
using System.Runtime.InteropServices;

public static class QuickDataGravity
{
    [DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
    public static extern void generate_bodies(IntPtr ctx, int count);

    [DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
    public static extern void get_body_positions(IntPtr ctx, float[] target, int count);

    [DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
    public static extern void get_body_radii(IntPtr ctx, float[] target, int count);

	[DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
	public static extern void set_sim_timestep(IntPtr ctx, float timestep);

    [DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
    public static extern void simulate(IntPtr ctx);

    [DllImport("QuickData.Gravity", CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr start_new();
}
