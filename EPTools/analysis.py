from .utils import *
import sys, os, subprocess
from pathlib import Path

__all__ = ['extract_curve_sh','get_clip_stamp','subtract_curve','make_specs_sh','make_rsp_sh']

def extract_curve_sh(evtfile, src_reg, bkg_reg, binsize=1, root='./', out_dir='./', output_name='xselect_lc'):
    
    if not Path(out_dir).is_absolute():
        out_dir = os.path.join(root,out_dir)
    if not os.path.exists(out_dir):
        print(f"Creating path {out_dir}")
        os.mkdir(out_dir)
    script = f"""
#!/bin/bash
rm -rf {os.path.join(out_dir,"src.lc")}
rm -rf {os.path.join(out_dir,"bkg.lc")}

# ================================
# User settings
# ================================
evtfile="{evtfile}"      # event file
binsize={binsize}              # time bin size (seconds)

# ================================
# Create xselect script
# ================================
cat <<EOF > xsel_lc.xco
xsel
read events
./
$evtfile
yes

# --------- basic filtering ----------
filter grade 0-12
select event "status==b0"

# --------- SOURCE light curve ----------
filter region {src_reg}
extract curve binsize=$binsize
save curve {os.path.join(out_dir,"src.lc")}
clear region
filter region {bkg_reg}
extract curve binsize=$binsize
save curve {os.path.join(out_dir,"bkg.lc")}

exit
no
EOF

# ================================
# Run xselect
# ================================
xselect < xsel_lc.xco

echo "Light curves saved in $outdir/"
    """
    if not output_name.endswith('.sh'):
        output_name = output_name+".sh"
    sh_path = os.path.join(root,output_name)
    with open(sh_path, "w") as f:
        f.write(script)

    os.chmod(sh_path, 0o755)
    print(f"xselect_lc.sh is stored in {sh_path}")
    
    
def get_clip_stamp(net_lc, ncounts, trange=None, out_file=None):
    """
    Return time bins [t_start, t_end] such that each bin contains ~ncounts.

    Parameters
    ----------
    net_lc : str
        FITS light curve file
    ncounts : float
        Counts per segment

    Returns
    -------
    segments : np.ndarray, shape (N, 2)
        [[t_start, t_end], ...]
    """


    with fits.open(net_lc) as hdul:

        # ---- Header ----
        hdr = hdul[0].header
        tstart = hdr.get('TSTART', 0.0)
        tend = hdr.get("TSTOP",0.0)

        # ---- Data ----
        data = hdul[1].data
        time = data['TIME']
        rate = data['RATE']
        
        # fliter time within trange
        if trange is None:
            trange = [0, tend-tstart]
        else:
            idx = (time > trange[0]) & (time < trange[1])
            time = time[idx]
            rate = rate[idx]
            
        

        # ---- bin width ----
        if 'TIMEDEL' in data.columns.names:
            dt = data['TIMEDEL']
        else:
            dt = np.median(np.diff(time))

        counts = rate * dt

        # ---- 累积 counts ----
        cum_counts = np.cumsum(counts)

        segments = []
        next_threshold = ncounts

        # 当前段起点
        t0 = time[0]

        for i in range(len(cum_counts)):
            if cum_counts[i] >= next_threshold:

                t1 = time[i]

                segments.append([
                    t0 + tstart,
                    t1 + tstart
                ])

                # 更新下一段
                t0 = time[i]
                next_threshold += ncounts
        if segments[-1][1] < trange[1]:
            segments.append([segments[-1][1],tend])

        if out_file is None:
            return np.array(segments)
        else:
            np.savez(out_file,np.array(segments))
    
    

    
def subtract_curve(src_lc, bkg_lc, alpha, out_file=None):
    """
    Subtract background light curve from source light curve.

    net = src - alpha * bkg

    Parameters
    ----------
    src_lc : str
        Source light curve FITS file
    bkg_lc : str
        Background light curve FITS file
    alpha : float
        Scaling factor (e.g., area ratio)
    out_file : str, optional
        Output FITS file

    Returns
    -------
    net_hdul : fits.HDUList
    """

    # -------------------------
    # Read files
    # -------------------------
    with fits.open(src_lc) as src_hdul, fits.open(bkg_lc) as bkg_hdul:

        src_data = src_hdul[1].data
        bkg_data = bkg_hdul[1].data

        # -------------------------
        # Check TIME alignment
        # -------------------------
        if not np.allclose(src_data['TIME'], bkg_data['TIME'], atol=1e-6):
            raise ValueError("Source and background TIME arrays do not match!")

        time = src_data['TIME']

        # -------------------------
        # RATE subtraction
        # -------------------------
        src_rate = src_data['RATE']
        bkg_rate = bkg_data['RATE']

        net_rate = src_rate - alpha * bkg_rate

        # -------------------------
        # ERROR propagation
        # -------------------------
        if 'ERROR' in src_data.columns.names and 'ERROR' in bkg_data.columns.names:
            src_err = src_data['ERROR']
            bkg_err = bkg_data['ERROR']

            net_err = np.sqrt(src_err**2 + (alpha * bkg_err)**2)
        else:
            net_err = None

        # -------------------------
        # Create output HDU
        # -------------------------
        net_hdul = src_hdul.copy()

        net_data = net_hdul[1].data

        net_data['RATE'] = net_rate

        if net_err is not None:
            if 'ERROR' in net_data.columns.names:
                net_data['ERROR'] = net_err

        # -------------------------
        # Optional: update header
        # -------------------------
        net_hdul[1].header['HISTORY'] = 'Background subtracted'
        net_hdul[1].header['HISTORY'] = f'alpha = {alpha}'

        # -------------------------
        # Save file
        # -------------------------
        if out_file is not None:
            net_hdul.writeto(out_file, overwrite=True)
            print(f"[OK] Saved to {out_file}")

    return net_hdul
    
    
def make_specs_sh(
    evtfile,
    stamp_file,
    root='./',
    out_dir='./',
    src_reg='src.reg',
    bkg_reg='bkg.reg',
    script_name='make_specs.sh'
):
    """
    Generate a bash script to extract time-resolved spectra using xselect.

    Parameters
    ----------
    evtfile : str
        Event file name (relative to root)
    stamp_file : str
        Numpy file storing [[t0, t1], ...]
    root : str
        Directory containing evt and region files
    out_dir : str
        Directory to save script
    src_reg : str
        Source region file
    bkg_reg : str
        Background region file
    script_name : str
        Output shell script name
    """

    # -------------------------
    # Load time stamps
    # -------------------------
    time_stamp = np.load(stamp_file)['arr_0']  # shape (N, 2)

    # -------------------------
    # Path handling
    # -------------------------
    root = Path(root).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    clipped_dir = out_dir / "int00"
    
    if Path(evtfile).is_absolute():
        evtfile = evtfile.split('/')[-1]

    rel_dir = root.relative_to(clipped_dir) if root.is_relative_to(clipped_dir) else Path(
        os.path.relpath(root, clipped_dir)
    )

    # -------------------------
    # Build time_ranges string
    # -------------------------
    time_ranges_str = ""
    for item in time_stamp:
        t0, t1 = item[0], item[1]
        time_ranges_str += f'"{t0}, {t1}"\n'


    # -------------------------
    # Build script
    # -------------------------
    script = f"""#!/bin/bash

time_ranges=(
{time_ranges_str}
)

for i in "${{!time_ranges[@]}}"; do
    folder=$(printf "int%02d" "$i")
    mkdir -p "$folder"
    cd "$folder" || exit

    cat <<EOF > xsel_script.xco
xsel
read events
{rel_dir}/
{evtfile}
yes
filter grade 0-12
select event "status==b0"
filter time scc
${{time_ranges[i]}}
x
filter region {rel_dir}/{src_reg}
extract spectrum
save spectrum src.pha
extract events
save events src.fits
no
clear events
clear region
filter region {rel_dir}/{bkg_reg}
extract spectrum
save spectrum bkg.pha
extract events
save events bkg.fits
no
exit
no
EOF

    xselect < xsel_script.xco

    cd ..
done
"""

    # -------------------------
    # Write file
    # -------------------------
    script_path = out_dir / script_name
    with open(script_path, 'w') as f:
        f.write(script)

    # 赋执行权限
    script_path.chmod(0o755)

    print(f"[OK] Script saved to {script_path}")

    return script_path


def make_rsp_sh(
    stamp_file,
    root='./',
    out_dir='./',
    expfile='expo.fits',
    prefix='src',
    script_name='make_rsp.sh'
):
    """
    Generate a bash script to create ARF and RMF for each time segment.

    Parameters
    ----------
    stamp_file : str
        Numpy file storing [[t0, t1], ...]
    root : str
        Directory containing exposure file
    out_dir : str
        Directory where intXX folders exist
    expfile : str
        Exposure map file name (relative to root)
    prefix : str
        Prefix of spectrum files (e.g., fxt_a → fxt_a_src.pha)
    script_name : str
        Output script name
    """

    # -------------------------
    # Load time stamps
    # -------------------------
    time_stamp = np.load(stamp_file)['arr_0']

    # -------------------------
    # Path handling
    # -------------------------
    root = Path(root).resolve()
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    clipped_dir = out_dir / "int00"

    # 相对路径（从 intXX 指向 root）
    rel_dir = root.relative_to(clipped_dir) if root.is_relative_to(clipped_dir) else Path(
        os.path.relpath(root, clipped_dir)
    )

    # -------------------------
    # Build script
    # -------------------------
    time_ranges_str = ""
    for item in time_stamp:
        t0, t1 = item[0], item[1]
        time_ranges_str += f'"{t0}, {t1}"\n'


    # -------------------------
    # Build script
    # -------------------------
    script = f"""#!/bin/bash

time_ranges=(
{time_ranges_str}
)
    
for i in "${{!time_ranges[@]}}"; do
    folder=$(printf "int%02d" "$i")
    echo "Processing $folder..."

    cd "$folder" || {{ echo "Failed to enter $folder"; exit 1; }}

    fxtarfgen specfile={prefix}.pha expfile={rel_dir}/{expfile} outfile={prefix}.arf
    fxtrmfgen specfile={prefix}.pha outfile={prefix}.rmf

    cd ..
done
"""

    # -------------------------
    # Write script
    # -------------------------
    script_path = out_dir / script_name
    with open(script_path, 'w') as f:
        f.write(script)

    script_path.chmod(0o755)

    print(f"[OK] Script saved to {script_path}")

    return script_path