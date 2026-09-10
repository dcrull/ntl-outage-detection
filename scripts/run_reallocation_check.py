#!/usr/bin/env python3
"""Reproduce original built-form allocation and test native-cell equivalence."""
import argparse
from dataclasses import asdict
from datetime import date
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from rasterio.transform import from_origin
from rasterio.warp import transform_bounds
import yaml
from scipy import sparse
from pyproj import Transformer

ROOT=Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path: sys.path.insert(0,str(ROOT))
from scripts.ingest_viirs import GROUP, cog, digest, now, save_json
from scripts.run_counterfactual import load_cache, evaluate
from scripts.counterfactual import Settings, fit_models, forecast, evidence
from scripts.reallocation import upstream_modules, native_overlap_matrix, fine_source_indices, aggregate_complete, allocate_nightly


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--upstream',type=Path,default=Path.home()/'research/ntl-psf-disaggregation')
    parser.add_argument('--work',type=Path,default=ROOT/'data/reallocation')
    parser.add_argument('--output',type=Path,default=ROOT/'data/published/reallocation_check')
    parser.add_argument('--prepare-only',action='store_true')
    args=parser.parse_args()
    work,out=args.work,args.output
    work.mkdir(parents=True,exist_ok=True); out.mkdir(parents=True,exist_ok=True)
    source_root=ROOT/'data/viirs/VNP46A2.002'
    data=load_cache(source_root,date(2025,4,2))
    cells=data['cells']; county=data['county']; days=data['dates']; cut=data['cutoff']
    modules=upstream_modules(args.upstream)
    op,kernels,bundle=(modules[k] for k in ('operator','kernels','overture_bundle'))
    config=yaml.safe_load((args.upstream/'configs/psf_disaggregation.yaml').read_text())
    built=config['allocation_proxies']['built_form']; release=built['overture_release']
    source_paths=[work/'overture'/release/f'{name}.parquet' for name in ('buildings','segments')]
    for path in source_paths:
        if not path.exists(): raise FileNotFoundError(f'Retrieve pinned Overture input: {path}')
        if pq.ParquetFile(path).metadata.num_rows<1: raise ValueError(f'Empty source {path}')
    reference=data['transform']; h,w=county.shape
    bounds=(reference.c,reference.f+h*reference.e,reference.c+w*reference.a,reference.f)
    b=transform_bounds('EPSG:4326','EPSG:32616',*bounds,densify_pts=101)
    b=tuple((math.floor(v/10)*10-1100 if i<2 else math.ceil(v/10)*10+1100) for i,v in enumerate(b))
    shape=(round((b[3]-b[1])/10),round((b[2]-b[0])/10)); transform=from_origin(b[0],b[3],10,10)
    bbox=transform_bounds('EPSG:32616','EPSG:4326',*b,densify_pts=101)
    upstream_files=list((args.upstream/'src/nocturne/disaggregate').glob('*.py'))
    upstream_files += [args.upstream/'configs/psf_disaggregation.yaml',args.upstream/'docs/psf-disaggregation-results.md']
    contract={'upstream_commit':subprocess.check_output(['git','-C',str(args.upstream),'rev-parse','HEAD'],text=True).strip(),
              'upstream_sha256':{str(p.relative_to(args.upstream)):digest(p) for p in upstream_files},
              'structural_inputs':{str(p.resolve()):digest(p) for p in source_paths},
              'resolution_m':10,'crs':'EPSG:32616','shape':list(shape),'bounds':list(b), 'bbox':list(bbox),
              'overture_release':release,'proxy':'built_form_primary','water_prior':'none',
              'kernel_radius_m':500,'proxy_floor':.05,'denominator_epsilon_relative':1e-6,
              'minimum_kernel_support':1.,'minimum_native_footprint_support':1.,
              'support_tolerance':1e-6,'target_cells_sha256':digest(ROOT/'data/published/counterfactual_dyer/cells.csv'),
              'temporal_caveat':'Contemporary 2026 structural proxy for 2025 radiance, matching upstream release; not historical reconstruction.'}
    # The CLI sidecars record the actual release and spatial extent requested.
    for p in source_paths:
        state_path=Path(str(p)+'.state')
        if not state_path.exists(): raise FileNotFoundError(f'Missing Overture provenance: {state_path}')
        state=json.loads(state_path.read_text())
        if state['last_release']!=release: raise ValueError('Overture release mismatch')
        sb=state['bbox']; expected=dict(zip(('xmin','ymin','xmax','ymax'),bbox))
        for key in expected:
            if (key in ('xmin','ymin') and sb[key]>expected[key]+1e-8) or (key in ('xmax','ymax') and sb[key]<expected[key]-1e-8):
                raise ValueError('Overture download does not cover processing halo')
        contract['structural_inputs'][str(state_path.resolve())]=digest(state_path)
    static_path=work/'static_manifest.json'
    cached=False
    if static_path.exists():
        old=json.loads(static_path.read_text())
        if old['contract']!=contract: raise ValueError('Static contract changed; use separate --work directory')
        cached=all((work/name).exists() and digest(work/name)==sha for name,sha in old['outputs'].items())
    # Derive the raw crop from the entire projected workspace, including its
    # halo and rotated corners, rather than the county-center bounding box.
    row0=math.floor((40-bbox[3])*240)-1
    col0=math.floor((bbox[0]+90)*240)-1
    row1=math.ceil((40-bbox[1])*240)+1
    col1=math.ceil((bbox[2]+90)*240)+1
    if min(row0,col0)<0 or max(row1,col1)>2400:
        raise ValueError('Projected workspace exceeds the cached h09v05 tile')
    source_shape=(row1-row0,col1-col0)
    kernel=kernels.circular_mean_kernel(radius_m=500,resolution_m=10)
    corr=op._correlate_constant_zero
    if not cached:
        print(f'Preparing {shape} projected workspace, original 10 m structural formula...',flush=True)
        tr=Transformer.from_crs('EPSG:4326','EPSG:32616',always_xy=True)
        bf,bmetrics=bundle._rasterize_building_fraction(source_paths[0],shape=shape,transform=transform,
                    source_to_target=tr,target_resolution_m=10,subpixel_resolution_m=2,progress_label='Dyer')
        road,rmetrics=bundle._rasterize_weighted_road_length(source_paths[1],shape=shape,transform=transform,
                    source_to_target=tr,road_weights=built['road_class_weights'],unlisted_road_class_policy='error',
                    maximum_segment_length_m=1,conservation_relative_tolerance=.0001,progress_label='Dyer')
        density=np.clip(road.astype('float64')/100*1e6/20000,0,1).astype('float32')
        raw=(.7*np.sqrt(bf)+.3*density).astype('float32')
        proxy=(.05+.95*raw).astype('float32')
        del bf,road,density,raw
        normalized,normalization=op.normalize_proxy_mean_one(proxy)
        denominator=corr(normalized,kernel.weights)
        geom=corr(np.ones(shape,dtype='float32'),kernel.weights)
        denominator=np.divide(denominator,geom,out=np.full(shape,np.nan,dtype='float32'),where=geom>1e-7)
        gain=op.allocation_gain_from_components(normalized,denominator,valid_mask=geom>=1-1e-6,
                                               denominator_epsilon=1e-6*normalization.mean_after)
        del proxy,normalized,denominator
        np.save(work/'gain.npy',gain); np.save(work/'geometric_support.npy',geom)
        del gain,geom
        mapping=fine_source_indices(transform,shape,row0,col0,source_shape)
        np.save(work/'fine_source_index.npy',mapping); del mapping
        overlap,areas=native_overlap_matrix(cells,transform,shape,progress=True)
        sparse.save_npz(work/'native_overlap.npz',overlap); np.save(work/'native_area_m2.npy',areas)
        del overlap
        names=['gain.npy','geometric_support.npy','fine_source_index.npy','native_overlap.npz','native_area_m2.npy']
        save_json(static_path,{'contract':contract,'normalization':normalization.to_metadata(),
                  'building_rasterization':bmetrics,'road_rasterization':rmetrics,
                  'outputs':{name:digest(work/name) for name in names}})
    if args.prepare_only: print('Static workspace ready',flush=True); return
    gain=np.load(work/'gain.npy',mmap_mode='r'); geom=np.load(work/'geometric_support.npy',mmap_mode='r')
    mapping=np.load(work/'fine_source_index.npy',mmap_mode='r'); W=sparse.load_npz(work/'native_overlap.npz')
    areas=np.load(work/'native_area_m2.npy')
    stage={'status':'running','started_utc':now(),'contract':contract,'static_sha256':digest(static_path),
           'input_sha256':data['hashes'],'code_sha256':{str(p.relative_to(ROOT)):digest(p) for p in
             [Path(__file__),ROOT/'scripts/reallocation.py',ROOT/'scripts/counterfactual.py',ROOT/'scripts/run_counterfactual.py']}}
    save_json(out/'manifest.json',stage)
    nightly=work/'nightly'; nightly.mkdir(exist_ok=True)
    records=[]; full={'direct':data['arrays']['observed']}; variants=['upsample','uniform','allocated']
    for key in variants+['allocation_support']: full[key]=np.full((len(days),len(cells)),np.nan)
    for i,day in enumerate(days):
        begin=time.monotonic(); cache=nightly/f'{day}.npz'; receipt=nightly/f'{day}.json'
        src=data['sources'][i]; raw_sha=None
        if src['source_available']:
            raw_path=source_root/'raw'/src['source_filename']; raw_sha=digest(raw_path)
            source_receipt=json.loads(Path(str(raw_path)+'.json').read_text())
            if raw_sha!=source_receipt['sha256']: raise ValueError('Raw halo source hash mismatch')
        key={'static_sha256':stage['static_sha256'],'source_sha256':raw_sha,
             'code_sha256':stage['code_sha256'],'date':day,'strict_QA':'same_as_direct_A2'}
        reuse=receipt.exists() and cache.exists()
        if reuse:
            previous=json.loads(receipt.read_text()); reuse=previous['key']==key and previous['sha256']==digest(cache)
        if reuse:
            with np.load(cache) as values: result={name:values[name] for name in variants+['allocation_support']}
        else:
            if src['source_available']:
                with h5py.File(raw_path,'r') as f:
                    g=f[GROUP]; sl=(slice(row0,row0+source_shape[0]),slice(col0,col0+source_shape[1]))
                    rad=g['DNB_BRDF-Corrected_NTL'][sl].astype('float32')
                    fill=float(np.asarray(g['DNB_BRDF-Corrected_NTL'].attrs['_FillValue']).ravel()[0])
                    cm=g['QF_Cloud_Mask'][sl]; mq=g['Mandatory_Quality_Flag'][sl]; snow=g['Snow_Flag'][sl]
                    good=np.isfinite(rad)&(rad!=fill)&(mq==0)&(cm!=65535)&(((cm>>6)&3)==0)&(((cm>>4)&3)==3)&((cm&1)==0)&(((cm>>9)&1)==0)&(snow==0)
                    rad=np.where(good,rad,np.nan)
                # Verify raw halo extraction against the existing direct input on every target cell.
                np.testing.assert_equal(rad[cells.source_row.to_numpy()-row0,cells.source_col.to_numpy()-col0],full['direct'][i].astype('float32'))
                fine=rad.ravel()[mapping]
            else: fine=np.full(shape,np.nan,dtype='float32')
            allocated,uniform=allocate_nightly(fine,gain,kernel,corr,geom)
            a,support=aggregate_complete(W,allocated)
            u,_=aggregate_complete(W,uniform); direct,_=aggregate_complete(W,fine)
            result={'allocated':a,'uniform':u,'upsample':direct,'allocation_support':support}
            np.savez_compressed(cache,**result); save_json(receipt,{'key':key,'sha256':digest(cache)})
            del fine,allocated,uniform
        for name,value in result.items(): full[name][i]=value
        valid=np.isfinite(full['direct'][i])&np.isfinite(result['allocated'])&np.isfinite(result['upsample'])&np.isfinite(result['uniform'])
        for variant in variants:
            error=result[variant][valid]-full['direct'][i,valid]
            all_direct=np.isfinite(full['direct'][i])
            rec={'date':day,'variant':variant,'direct_usable_cells':int(all_direct.sum()),'common_cells':int(valid.sum()),
                 'max_abs_error':float(np.max(np.abs(error))) if len(error) else None,
                 'mae':float(np.mean(np.abs(error))) if len(error) else None,
                 'rmse':float(np.sqrt(np.mean(error**2))) if len(error) else None,
                 'bias':float(np.mean(error)) if len(error) else None,
                 'area_weighted_bias':float(np.average(error,weights=areas[valid])) if len(error) else None,
                 'fraction_numerically_equal':float(np.mean(np.isclose(result[variant][valid],full['direct'][i,valid],rtol=1e-5,atol=1e-5))) if len(error) else None}
            records.append(rec)
        print(f'{day}: common {valid.sum()}/{np.isfinite(full["direct"][i]).sum()}, {time.monotonic()-begin:.1f}s'+(' cached' if reuse else ''),flush=True)
    # Identical cell-night masks for downstream comparisons; cloud/support loss is explicit.
    common=np.logical_and.reduce([np.isfinite(full[k]) for k in ['direct']+variants])
    settings=Settings(); model_arrays={}; validations=[]; model_differences=[]
    outfiles=[]
    def table(name,df): df.to_csv(out/name,index=False,float_format='%.9g'); outfiles.append(name)
    def raster(name,values,units='nW cm-2 sr-1'):
        grid=np.full(county.shape,np.nan,dtype='float32'); grid[county]=values
        cog(out/name,grid,reference,np.nan,units=units,description='Native reallocation comparison; not outage classification')
        outfiles.append(name)
    table('native_equivalence_by_date.csv',pd.DataFrame(records)); table('cells.csv',cells.assign(area_m2=areas))
    for variant in ['direct']+variants:
        matched=np.where(common,full[variant],np.nan)
        fitted=fit_models(matched[:cut],settings)
        summary,_,_=evaluate(matched,cut,settings,days); summary['input_variant']=variant
        validations.append(summary)
        model_arrays[f'{variant}_training_count']=fitted['count']
        for method in ('M0','M1'):
            pred=forecast(fitted,method,np.arange(1,len(days)-cut+1)); scores=evidence(pred,matched[cut:])
            for name,v in {**pred,**scores}.items(): model_arrays[f'{variant}_{method}_{name}']=v
            for j,day in enumerate(days[cut:]):
                for name in ('center','q05','q95','residual','lower_tail'):
                    values=pred[name][j] if name in pred else scores[name][j]
                    raster(f'{variant}/{method}/{day}/{name}.tif',values,'dimensionless' if name=='lower_tail' else 'nW cm-2 sr-1')
            if variant!='direct':
                for name,v in {**pred,**scores}.items():
                    ref=model_arrays[f'direct_{method}_{name}']; ok=np.isfinite(v)&np.isfinite(ref); e=v[ok]-ref[ok]
                    model_differences.append({'variant':variant,'method':method,'field':name,'n':len(e),
                            'mae':float(np.mean(np.abs(e))) if len(e) else None,'max_abs_error':float(np.max(np.abs(e))) if len(e) else None})
        for i,day in enumerate(days):
            raster(f'native/{variant}/{day}.tif',full[variant][i])
    table('validation_summary.csv',pd.concat(validations)); table('model_differences.csv',pd.DataFrame(model_differences))
    np.savez_compressed(out/'comparison_arrays.npz',**full,**model_arrays,common=common,dates=np.array(days),cutoff=cut,county=county,transform=np.array(list(reference)[:6]));outfiles.append('comparison_arrays.npz')
    table('support_by_date.csv',pd.DataFrame({'date':days,'direct_usable':np.isfinite(full['direct']).sum(axis=1),'common_usable':common.sum(axis=1)}))
    stage.update(status='complete',finished_utc=now(),settings=asdict(settings),
            interpretation='Tests of equality, not enforced conservation. Strict nightly QA retained; complete fine-footprint support required. No TING.',
            output_sha256={name:digest(out/name) for name in outfiles})
    save_json(out/'manifest.json',stage)
    print('Complete:',out,flush=True)

if __name__=='__main__': main()
