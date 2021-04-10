import os
import sys
from lsaBGC.classes.Pan import Pan
class GCF(Pan):
  def __init__(self, bgc_genbanks_listing, gcf_id='GCF', lineage_name='Unnamed lineage'):
    super().__init__(bgc_genbanks_listing, lineage_name=lineage_name)
    self.gcf_id = gcf_id

