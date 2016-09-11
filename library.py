from qgis.core import *
from PyQt4.QtCore import QVariant
import math

#from eq_distant_dialog import EqDistantDialog


class Library:
    def __init__(self, layer_a, layer_b, intv, jarak_klaim=0):
        self.layer_a = layer_a
        self.layer_b = layer_b
        self.intv = intv
        self.crs = layer_a.crs().authid()
        self.jarak_klaim = jarak_klaim

    def konversi_garis_ke_titik(self, list_ft_garis, attr):
        # create_poi
        layer_ttk = QgsVectorLayer("Point?crs=" + self.crs, "Point Layer", "memory")
        provider_layer_ttk = layer_ttk.dataProvider()
        provider_layer_ttk.addAttributes([QgsField("id", QVariant.Int), QgsField("ket", QVariant.String)])
        list_ft_ttk = []
        fid = 0
        for ft_garis in list_ft_garis:
            list_koord_garis = ft_garis.geometry().asPolyline()
            for koord in list_koord_garis:
                ttk = QgsPoint(koord[0], koord[1])
                geom_ttk = QgsGeometry.fromPoint(ttk)
                feat_ttk = QgsFeature()
                feat_ttk.setGeometry(geom_ttk)
                feat_ttk.setAttributes([fid, attr])
                fid += 1
                list_ft_ttk.append(feat_ttk)
        provider_layer_ttk.addFeatures(list_ft_ttk)
        layer_ttk.startEditing()
        layer_ttk.commitChanges()
        return layer_ttk

    def cari_titik_tengah(self, ttk_a, ttk_b):
        xa = ttk_a.x()
        ya = ttk_a.y()
        xb = ttk_b.x()
        yb = ttk_b.y()
        xt = xa - (xa-xb)/2
        yt = ya - (ya-yb)/2
        ttk_tengah = QgsPoint(xt, yt)
        return ttk_tengah

    def hitung_circumcenter(self, ttk_a, ttk_b, ttk_c):
        xa = ttk_a.x()
        ya = ttk_a.y()
        xb = ttk_b.x()
        yb = ttk_b.y()
        xc = ttk_c.x()
        yc = ttk_c.y()
        a2 = pow(xa, 2) + pow(ya, 2)
        b2 = pow(xb, 2) + pow(yb, 2)
        c2 = pow(xc, 2) + pow(yc, 2)
        # ------------------
        konstanta = 2*(xa*(yb - yc) + xb*(yc - ya) + xc*(ya - yb))
        xo = (a2*(yb - yc) + b2*(yc - ya) + c2*(ya - yb))/konstanta
        yo = (a2*(xc - xb) + b2*(xa - xc) + c2*(xb - xa))/konstanta
        ttk_o = QgsPoint(xo, yo)
        return ttk_o

    def cari_jarak_titik(self, ttk_a, ttk_b):
        return math.sqrt(ttk_a.sqrDist(ttk_b))

    def cari_sisi_garis(self, ttk_a, ttk_b, ttk_sisi):
        xa = ttk_a.x()
        ya = ttk_a.y()
        xb = ttk_b.x()
        yb = ttk_b.y()
        x_sisi = ttk_sisi.x()
        y_sisi = ttk_sisi.y()
        s = (x_sisi - xa)*(yb - ya) - (y_sisi - ya)*(xb - xa)
        if s < 0:
            sisi = -1
        elif s > 0:
            sisi = 1
        else:
            sisi = 0
        return sisi

    def titik_pada_garis(self, titik, list_geom_garis):
        hasil = {}
        for geom_garis in list_geom_garis:
            jarak, t, int = geom_garis.closestSegmentWithContext(titik)
            hasil[t] = jarak
        ttk_terdekat = min(hasil, key=hasil.get)
        return ttk_terdekat

    # Opposite
    def buat_garis_tgk_lurus(self, jarak, ttk_a, ttk_b, ttk_arah_garis):
        ttk_tgh = self.cari_titik_tengah(ttk_a, ttk_b)
        xt = ttk_tgh.x()
        yt = ttk_tgh.y()
        azm = ttk_a.azimuth(ttk_b)
        sisi = self.cari_sisi_garis(ttk_a, ttk_b, ttk_arah_garis)
        x_tgk = 0
        y_tgk = 0
        if azm > 0 and azm <= 90:
            if sisi == 1:
                x_tgk = xt + jarak*math.cos(math.radians(azm))
                y_tgk = yt - jarak*math.sin(math.radians(azm))
            elif sisi == -1:
                x_tgk = xt - jarak*math.cos(math.radians(azm))
                y_tgk = yt + jarak*math.sin(math.radians(azm))
        elif azm > 90 and azm <= 180:
            if sisi == 1:
                x_tgk = xt + jarak*math.cos(math.radians(azm))
                y_tgk = yt - jarak*math.sin(math.radians(azm))
            elif sisi == -1:
                x_tgk = xt - jarak*math.cos(math.radians(azm))
                y_tgk = yt + jarak*math.sin(math.radians(azm))
        elif azm < 0 and azm >= -90:
            if sisi == 1:
                x_tgk = xt + jarak*math.cos(math.radians(azm))
                y_tgk = yt - jarak*math.sin(math.radians(azm))
            elif sisi == -1:
                x_tgk = xt - jarak*math.cos(math.radians(azm))
                y_tgk = yt + jarak*math.sin(math.radians(azm))
        elif azm < -90 and azm >= -180:
            if sisi == 1:
                x_tgk = xt + jarak*math.cos(math.radians(azm))
                y_tgk = yt - jarak*math.sin(math.radians(azm))
            elif sisi == -1:
                x_tgk = xt - jarak*math.cos(math.radians(azm))
                y_tgk = yt + jarak*math.sin(math.radians(azm))
        else:
            raise ValueError('Azimuth Error')
        ttk_tgk = QgsPoint(x_tgk, y_tgk)
        garis_tgk = QgsGeometry.fromPolyline([ttk_tgh, ttk_tgk])
        return garis_tgk

    def cari_fitur_terdekat(self, ttk, list_ft):
        jarak_terdekat = 999999999
        fitur_terdekat = QgsFeature()
        for ft in list_ft:
            ttk_ft = ft.geometry().asPoint()
            jarak = math.sqrt(ttk.sqrDist(ttk_ft))
            if jarak < jarak_terdekat:
                jarak_terdekat = jarak
                fitur_terdekat = ft
        return fitur_terdekat

    def filter_titik(self, ttk_a, ttk_b, ttk_akhir_a, ttk_akhir_b, list_ft):
        list_baru = []
        ttk_tgh_aw = self.cari_titik_tengah(ttk_a, ttk_b)
        ttk_tgh_ak = self.cari_titik_tengah(ttk_akhir_a, ttk_akhir_b)
        sisi_aw = self.cari_sisi_garis(ttk_a, ttk_b, ttk_tgh_ak)
        sisi_ak = self.cari_sisi_garis(ttk_akhir_a, ttk_akhir_b, ttk_tgh_aw)
        for ft in list_ft:
            ttk_ft = ft.geometry().asPoint()
            sisi_ft_aw = self.cari_sisi_garis(ttk_a, ttk_b, ttk_ft)
            sisi_ft_ak = self.cari_sisi_garis(ttk_akhir_a, ttk_akhir_b, ttk_ft)
            if sisi_ft_aw == sisi_aw and sisi_ft_ak == sisi_ak:
                list_baru.append(ft)
        return list_baru

    def cek_perpotongan(self, ttk_a, ttk_b, jarak, ttk_eq, sisi_ak, list_ft):
        hasil = []
        #ft_berikutnya = QgsFeature()
        for ft in list_ft:
            ttk_ft = ft.geometry().asPoint()
            jarak_ttk_ft = math.sqrt(ttk_ft.sqrDist(ttk_eq))
            if jarak_ttk_ft < jarak and ttk_ft not in [ttk_a, ttk_b]:
                hasil.append(ft)
        if len(hasil) == 1:
            ft_berikutnya = hasil[0]
        elif len(hasil) > 1:
            ft_berikutnya = self.cari_fitur_terdekat(ttk_eq, hasil)
        else:
            ft_berikutnya = 0
        return ft_berikutnya

    def cari_fitur_ketiga(self, ttk_a, ttk_b, sisi_ak, list_ft, garis_tgk):
        jarak_interpolasi = 0
        while jarak_interpolasi <= garis_tgk.length():
            geom_eq = garis_tgk.interpolate(jarak_interpolasi)
            ttk_eq = geom_eq.asPoint()
            # cari jarak
            jarak_a = math.sqrt(ttk_eq.sqrDist(ttk_a))
            jarak_b = math.sqrt(ttk_eq.sqrDist(ttk_b))
            jarak_r = (jarak_a + jarak_b)/2
            ft_ketiga = self.cek_perpotongan(ttk_a, ttk_b, jarak_r, ttk_eq, sisi_ak, list_ft)
            if ft_ketiga != 0:
                ttk_ketiga = ft_ketiga.geometry().asPoint()
                sisi_ketiga = self.cari_sisi_garis(ttk_a, ttk_b, ttk_ketiga)
                if sisi_ketiga == sisi_ak:
                    break
            jarak_interpolasi += self.intv
        else:
            ft_ketiga = None
            geom_eq = None
        return ft_ketiga, geom_eq

    def layer_titik_konstruksi(self, featlist_ttk_konstruksi):
        '''fungsi ini untuk membuat fitur titik konstruksi '''
        lay_ttk_kons = QgsVectorLayer("Point?crs=" + self.crs, "Construction Point", "memory")
        prov_lay_ttk_kons = lay_ttk_kons.dataProvider()
        prov_lay_ttk_kons.addAttributes([QgsField("id",QVariant.String),QgsField("x",QVariant.Int), QgsField("y", QVariant.Int)])
        prov_lay_ttk_kons.addFeatures(featlist_ttk_konstruksi)
        lay_ttk_kons.startEditing()
        lay_ttk_kons.commitChanges()
        return lay_ttk_kons

    def layer_titik_circumcenter(self, featlist_ttk_circumcenter):
        lay_ttk_circumcenter = QgsVectorLayer("Point?crs=" + self.crs, "Equidistant Point", "memory")
        prov_lay_ttk_cc = lay_ttk_circumcenter.dataProvider()
        prov_lay_ttk_cc.addAttributes([QgsField("id",QVariant.String),QgsField("x",QVariant.Int), QgsField("y", QVariant.Int)])
        prov_lay_ttk_cc.addFeatures(featlist_ttk_circumcenter)
        lay_ttk_circumcenter.startEditing()
        lay_ttk_circumcenter.commitChanges()
        return lay_ttk_circumcenter

    def proses_hdp(self, ttk_awal_a, ttk_awal_b, ttk_akhir_a, ttk_akhir_b, list_ft):
        # /////////////// PERSIAPAN VARIABEL
        featlist_titik_circumcenter = []
        featlist_titik_konstruksi = []
        featlist_garis_konstruksi = []
        id_circumcenter = 1
        id_ttk_konstruksi_a = 1
        id_ttk_konstruksi_b = 1
        # ------
        feat_ttk_iter_a = QgsFeature()
        feat_ttk_iter_a.setGeometry(QgsGeometry.fromPoint(ttk_awal_a))
        feat_ttk_iter_a.setAttributes(["A"+str(id_ttk_konstruksi_a),
                                       ttk_awal_a.x(),
                                       ttk_awal_a.y()])
        featlist_titik_konstruksi.append(feat_ttk_iter_a)
        # ------
        feat_ttk_iter_b = QgsFeature()
        feat_ttk_iter_b.setGeometry(QgsGeometry.fromPoint(ttk_awal_b))
        feat_ttk_iter_b.setAttributes(["B"+str(id_ttk_konstruksi_b),
                                       ttk_awal_b.x(),
                                       ttk_awal_b.y()])
        featlist_titik_konstruksi.append(feat_ttk_iter_b)
        # /////////////// PERHITUNGAN VARIABEL
        ttk_tgh_aw = self.cari_titik_tengah(ttk_awal_a, ttk_awal_b)
        ttk_tgh_ak = self.cari_titik_tengah(ttk_akhir_a, ttk_akhir_b)
        jarak_tgh = math.sqrt(ttk_tgh_aw.sqrDist(ttk_tgh_ak))
        sisi_ak = self.cari_sisi_garis(ttk_awal_a, ttk_awal_b, ttk_tgh_ak)
        ttk_iter_a = feat_ttk_iter_a.geometry().asPoint()
        ttk_iter_b = feat_ttk_iter_b.geometry().asPoint()
        garis_tegak = self.buat_garis_tgk_lurus(jarak_tgh, ttk_iter_a, ttk_iter_b, ttk_tgh_ak)
        feat_ketiga, geom_eq = self.cari_fitur_ketiga(ttk_iter_a, ttk_iter_b, sisi_ak, list_ft, garis_tegak)
        while feat_ketiga is not None:
            ttk_ketiga = feat_ketiga.geometry().asPoint()
            ttk_circumcenter = self.hitung_circumcenter(ttk_iter_a, ttk_iter_b, ttk_ketiga)
            feat_ttk_circumcenter = QgsFeature()
            garis_kontrol_a = QgsFeature()
            garis_kontrol_b = QgsFeature()
            garis_kontrol_c = QgsFeature()
            feat_ttk_circumcenter.setGeometry(QgsGeometry().fromPoint(ttk_circumcenter))
            feat_ttk_circumcenter.setAttributes([str(id_circumcenter), ttk_circumcenter.x(), ttk_circumcenter.y()])
            featlist_titik_circumcenter.append(feat_ttk_circumcenter)

            g_garis_kontrol_a = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_iter_a])
            g_garis_kontrol_b = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_iter_b])
            garis_kontrol_a.setGeometry(g_garis_kontrol_a)
            garis_kontrol_a.setAttributes([str(id_circumcenter)+"a",
                                           g_garis_kontrol_a.length(),
                                           feat_ttk_iter_a[0],
                                           id_circumcenter])
            garis_kontrol_b.setGeometry(g_garis_kontrol_b)
            garis_kontrol_b.setAttributes([str(id_circumcenter)+"b",
                                           g_garis_kontrol_b.length(),
                                           feat_ttk_iter_b[0],
                                           id_circumcenter])
            g_garis_kontrol_c = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_ketiga])
            garis_kontrol_c.setGeometry(g_garis_kontrol_c)
            if feat_ketiga["ket"]=="A":
                feat_ttk_iter_a = feat_ketiga
                id_ttk_konstruksi_a += 1
                feat_ttk_iter_a.setAttributes(["A"+str(id_ttk_konstruksi_a), ttk_ketiga.x(), ttk_ketiga.y()])
                featlist_titik_konstruksi.append(feat_ttk_iter_a)
                garis_kontrol_c.setAttributes([str(id_circumcenter)+"c",
                                           g_garis_kontrol_c.length(),
                                           feat_ttk_iter_a[0],
                                           id_circumcenter])
                ttk_iter_a = feat_ttk_iter_a.geometry().asPoint()
            elif feat_ketiga["ket"]=="B":
                feat_ttk_iter_b = feat_ketiga
                id_ttk_konstruksi_b += 1
                feat_ttk_iter_b.setAttributes(["B"+str(id_ttk_konstruksi_b), ttk_ketiga.x(), ttk_ketiga.y()])
                featlist_titik_konstruksi.append(feat_ttk_iter_b)
                garis_kontrol_c.setAttributes([str(id_circumcenter)+"c",
                                           g_garis_kontrol_c.length(),
                                           feat_ttk_iter_b[0],
                                           id_circumcenter])
                ttk_iter_b = feat_ttk_iter_b.geometry().asPoint()
            else:
                raise ValueError("Attribute Error")
            featlist_garis_konstruksi.append(garis_kontrol_a)
            featlist_garis_konstruksi.append(garis_kontrol_b)
            featlist_garis_konstruksi.append(garis_kontrol_c)
            id_circumcenter += 1
            ttk_tgh_iter = self.cari_titik_tengah(ttk_iter_a, ttk_iter_b)
            jarak_tgh_iter = math.sqrt(ttk_tgh_iter.sqrDist(ttk_tgh_ak))
            garis_tegak_iter = self.buat_garis_tgk_lurus(jarak_tgh_iter, ttk_iter_a, ttk_iter_b, ttk_tgh_ak)
            feat_ketiga, geom_eq = self.cari_fitur_ketiga(ttk_iter_a, ttk_iter_b, sisi_ak, list_ft, garis_tegak_iter)
        else:
            print("Proses Selesai")
        return featlist_titik_circumcenter, featlist_titik_konstruksi, featlist_garis_konstruksi
    # Adjacent
    def buat_titik_tgk_lurus(self, jarak, ttk_a, ttk_b, sisi_garis):
        ttk_tgh = self.cari_titik_tengah(ttk_a, ttk_b)
        xt = ttk_tgh.x()
        yt = ttk_tgh.y()
        azm = ttk_a.azimuth(ttk_b)
        if sisi_garis > 0:
            x_tgk = xt + jarak*math.cos(math.radians(azm))
            y_tgk = yt - jarak*math.sin(math.radians(azm))
        elif sisi_garis < 0:
            x_tgk = xt - jarak*math.cos(math.radians(azm))
            y_tgk = yt + jarak*math.sin(math.radians(azm))
        else:
            raise ValueError('Side value is null')
        return QgsPoint(x_tgk, y_tgk)

    def proses_sb(self, ttk_awal_a, ttk_awal_b, ttk_akhir, list_ft):
        # /////////////// PERSIAPAN VARIABEL
        featlist_titik_circumcenter = []
        featlist_titik_konstruksi = []
        featlist_garis_konstruksi = []
        id_circumcenter = 1
        id_ttk_konstruksi_a = 1
        id_ttk_konstruksi_b = 1
        # ------
        feat_ttk_iter_a = QgsFeature()
        feat_ttk_iter_a.setGeometry(QgsGeometry.fromPoint(ttk_awal_a))
        feat_ttk_iter_a.setAttributes(["A"+str(id_ttk_konstruksi_a),
                                       ttk_awal_a.x(),
                                       ttk_awal_a.y()])
        featlist_titik_konstruksi.append(feat_ttk_iter_a)
        # ------
        feat_ttk_iter_b = QgsFeature()
        feat_ttk_iter_b.setGeometry(QgsGeometry.fromPoint(ttk_awal_b))
        feat_ttk_iter_b.setAttributes(["B"+str(id_ttk_konstruksi_b),
                                       ttk_awal_b.x(),
                                       ttk_awal_b.y()])
        featlist_titik_konstruksi.append(feat_ttk_iter_b)
        # ////////////////// PERHITUNGAN VARIABEL
        ttk_tgh = self.cari_titik_tengah(ttk_awal_a, ttk_awal_b)
        sisi_ak = self.cari_sisi_garis(ttk_awal_a, ttk_awal_b, ttk_akhir)
        jarak_tgh_aw = math.sqrt(ttk_tgh.sqrDist(ttk_akhir))
        jarak_tgh_o = self.jarak_klaim - jarak_tgh_aw
        ttk_iter_a = feat_ttk_iter_a.geometry().asPoint()
        ttk_iter_b = feat_ttk_iter_b.geometry().asPoint()
        ttk_tgk_ak = self.buat_titik_tgk_lurus(jarak_tgh_aw,
                                               ttk_iter_a,
                                               ttk_iter_b,
                                               sisi_ak)
        ttk_tgk_o = self.buat_titik_tgk_lurus(jarak_tgh_o,
                                              ttk_iter_a,
                                              ttk_iter_b,
                                              -sisi_ak)
        grs_tegak = QgsGeometry().fromPolyline([ttk_tgk_o, ttk_tgk_ak])
        feat_ketiga, geom_eq = self.cari_fitur_ketiga(ttk_iter_a,
                                                      ttk_iter_b,
                                                      sisi_ak,
                                                      list_ft,
                                                      grs_tegak)
        print "initialization complete"
        while feat_ketiga is not None:
            print "---------------"
            print "starting iteration"
            ttk_ketiga = feat_ketiga.geometry().asPoint()
            ttk_circumcenter = self.hitung_circumcenter(ttk_iter_a, ttk_iter_b, ttk_ketiga)
            feat_ttk_circumcenter = QgsFeature()
            garis_kontrol_a = QgsFeature()
            garis_kontrol_b = QgsFeature()
            garis_kontrol_c = QgsFeature()
            feat_ttk_circumcenter.setGeometry(QgsGeometry().fromPoint(ttk_circumcenter))
            feat_ttk_circumcenter.setAttributes([str(id_circumcenter), ttk_circumcenter.x(), ttk_circumcenter.y()])
            featlist_titik_circumcenter.append(feat_ttk_circumcenter)
            #
            g_garis_kontrol_a = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_iter_a])
            garis_kontrol_a.setGeometry(g_garis_kontrol_a)
            garis_kontrol_a.setAttributes([str(id_circumcenter)+"a",
                                           g_garis_kontrol_a.length(),
                                           feat_ttk_iter_a[0],
                                           id_circumcenter])
            g_garis_kontrol_b = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_iter_b])
            garis_kontrol_b.setGeometry(g_garis_kontrol_b)
            garis_kontrol_b.setAttributes([str(id_circumcenter)+"b",
                                           g_garis_kontrol_b.length(),
                                           feat_ttk_iter_b[0],
                                           id_circumcenter])
            g_garis_kontrol_c = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_ketiga])
            garis_kontrol_c.setGeometry(g_garis_kontrol_c)
            print "checking attribute"
            if feat_ketiga["ket"] == "A":
                feat_ttk_iter_a = feat_ketiga
                id_ttk_konstruksi_a += 1
                feat_ttk_iter_a.setAttributes(["A"+str(id_ttk_konstruksi_a), ttk_ketiga.x(), ttk_ketiga.y()])
                featlist_titik_konstruksi.append(feat_ttk_iter_a)
                garis_kontrol_c.setAttributes([str(id_circumcenter)+"c",
                                           g_garis_kontrol_c.length(),
                                           feat_ttk_iter_a[0],
                                           id_circumcenter])
                ttk_iter_a = feat_ttk_iter_a.geometry().asPoint()
            elif feat_ketiga["ket"]=="B":
                feat_ttk_iter_b = feat_ketiga
                id_ttk_konstruksi_b += 1
                feat_ttk_iter_b.setAttributes(["B"+str(id_ttk_konstruksi_b), ttk_ketiga.x(), ttk_ketiga.y()])
                featlist_titik_konstruksi.append(feat_ttk_iter_b)
                garis_kontrol_c.setAttributes([str(id_circumcenter)+"c",
                                           g_garis_kontrol_c.length(),
                                           feat_ttk_iter_b[0],
                                           id_circumcenter])
                ttk_iter_b = feat_ttk_iter_b.geometry().asPoint()
            else:
                raise ValueError("Attribute Error")
            print "attribute checking complete"
            featlist_garis_konstruksi.append(garis_kontrol_a)
            featlist_garis_konstruksi.append(garis_kontrol_b)
            featlist_garis_konstruksi.append(garis_kontrol_c)
            id_circumcenter += 1
            ttk_tgh_iter = self.cari_titik_tengah(ttk_iter_a, ttk_iter_b)
            jarak_tgh_iter = math.sqrt(ttk_tgh_iter.sqrDist(ttk_akhir))
            ttk_tgk_iter = self.buat_titik_tgk_lurus(jarak_tgh_iter,ttk_iter_a,ttk_iter_b,sisi_ak)
            grs_tegak = QgsGeometry.fromPolyline([ttk_circumcenter, ttk_tgk_iter])
            feat_ketiga, geom_eq = self.cari_fitur_ketiga(ttk_iter_a,
                                                          ttk_iter_b,
                                                          sisi_ak,
                                                          list_ft,
                                                          grs_tegak)
        else:
            print "Proses Selesai"
        feat_ttk_akhir = QgsFeature()
        feat_ttk_akhir.setGeometry(QgsGeometry().fromPoint(ttk_akhir))
        feat_ttk_akhir.setAttributes([str(id_circumcenter), ttk_circumcenter.x(), ttk_circumcenter.y()])
        featlist_titik_circumcenter.append(feat_ttk_akhir)
        return featlist_titik_circumcenter, featlist_titik_konstruksi, featlist_garis_konstruksi

    # Operasi Pembuatan Layer
    def layer_titik(self, nama_layer, list_feat):
        layer_titik = QgsVectorLayer("point?crs="+self.crs, nama_layer, "memory")
        provider_layer = layer_titik.dataProvider()
        provider_layer.addAttributes([QgsField("ID", QVariant.String),
                                      QgsField("X", QVariant.Int),
                                      QgsField("Y", QVariant.Int)
                                      ])
        provider_layer.addFeatures(list_feat)
        layer_titik.startEditing()
        layer_titik.commitChanges()
        QgsMapLayerRegistry.instance().addMapLayer(layer_titik)

    def layer_garis(self, nama_layer, list_feat):
        layer_garis_kons = QgsVectorLayer("Linestring?crs="+self.crs, nama_layer, "memory")
        provider_layer = layer_garis_kons.dataProvider()
        provider_layer.addAttributes([QgsField("ID", QVariant.String),
                                      QgsField("Length", QVariant.Int),
                                      QgsField("id_cons_pt", QVariant.String),
                                      QgsField("id_eq_pt", QVariant.String)
                                      ])
        provider_layer.addFeatures(list_feat)
        layer_garis_kons.startEditing()
        layer_garis_kons.commitChanges()
        QgsMapLayerRegistry.instance().addMapLayer(layer_garis_kons)

    def buat_layer_titik(self, list_geom_titik):
        layer_titik = QgsVectorLayer("point?crs="+self.crs, "Equidistant Point", "memory")
        provider_layer = layer_titik.dataProvider()
        list_ft = []
        for geom in list_geom_titik:
            feat_ttk = QgsFeature()
            feat_ttk.setGeometry(geom)
            list_ft.append(feat_ttk)
        provider_layer.addFeatures(list_ft)
        QgsMapLayerRegistry.instance().addMapLayer(layer_titik)

    def buat_layer_garis_eq(self, geom_garis):
        layer_garis = QgsVectorLayer("linestRing?crs="+self.crs, "Eq Line", "memory")
        provider_layer = layer_garis.dataProvider()
        feat_garis = QgsFeature()
        feat_garis.setGeometry(geom_garis)
        provider_layer.addFeatures([feat_garis])
        QgsMapLayerRegistry.instance().addMapLayer(layer_garis)

    def buat_layer_garis_k(self, list_geom_garis):
        layer_grs_kontrol = QgsVectorLayer("linestrIng?crs="+self.crs, "Construction Line", "memory")
        provider_layer = layer_grs_kontrol.dataProvider()
        list_feat_garis =[]
        for geom in list_geom_garis:
            feat_grs = QgsFeature()
            feat_grs.setGeometry(geom)
            list_feat_garis.append(feat_grs)
        provider_layer.addFeatures(list_feat_garis)
        QgsMapLayerRegistry.instance().addMapLayer(layer_grs_kontrol)

    def konversi_titik_ke_garis(self, list_geom_ttk):
        # konversi list geometri titik menjadi feature
        list_ft = []
        for geom in list_geom_ttk:
            feat = QgsFeature()
            feat.setGeometry(geom)
            list_ft.append(feat)
        list_hsl = []
        feat_aw = list_ft[0]
        ttk_aw = feat_aw.geometry().asPoint()
        list_hsl.append(ttk_aw)
        list_ft.remove(feat_aw)
        while len(list_ft) > 0:
            feat_n = self.cari_fitur_terdekat(ttk_aw, list_ft)
            ttk_n = feat_n.geometry().asPoint()
            list_hsl.append(ttk_n)
            list_ft.remove(feat_n)
            ttk_aw = ttk_n
        # konversi list titik ke layer garis
        layer_garis = QgsVectorLayer("lineString?crs=" + self.crs, "Equidistant Line", "memory")
        prov_layer_garis = layer_garis.dataProvider()
        feat_garis = QgsFeature()
        geom_garis = QgsGeometry.fromPolyline(list_hsl)
        feat_garis.setGeometry(geom_garis)
        prov_layer_garis.addFeatures([feat_garis])
        QgsMapLayerRegistry.instance().addMapLayer(layer_garis)

    def densify_line(self, geom_garis, dens_intv):
        list_pt_geom = []
        list_koord_titik = geom_garis.asPolyline()
        for koord in list_koord_titik:
            titik = QgsPoint(koord[0], koord[1])
            geom_titik = QgsGeometry.fromPoint(titik)
            list_pt_geom.append(geom_titik)
        intv_dist = dens_intv
        while intv_dist < geom_garis.length() :
            geom_ttk = geom_garis.interpolate(intv_dist)
            list_pt_geom.append(geom_ttk)
            intv_dist += dens_intv
        self.konversi_titik_ke_garis(list_pt_geom)

    def densify_line_geom(self, geom_garis, dens_intv):
        list_pt_geom = []
        list_koord_titik = geom_garis.asPolyline()
        for koord in list_koord_titik:
            titik = QgsPoint(koord[0], koord[1])
            geom_titik = QgsGeometry.fromPoint(titik)
            list_pt_geom.append(geom_titik)
        intv_dist = dens_intv
        while intv_dist < geom_garis.length() :
            geom_ttk = geom_garis.interpolate(intv_dist)
            list_pt_geom.append(geom_ttk)
            intv_dist += dens_intv
        list_ft = []
        for geom in list_pt_geom:
            feat = QgsFeature()
            feat.setGeometry(geom)
            list_ft.append(feat)
        feat_aw = list_ft[0]
        ttk_aw = feat_aw.geometry().asPoint()
        list_hsl = []
        list_hsl.append(ttk_aw)
        list_ft.remove(feat_aw)
        while len(list_ft) > 0:
            feat_n = self.cari_fitur_terdekat(ttk_aw, list_ft)
            ttk_n = feat_n.geometry().asPoint()
            list_hsl.append(ttk_n)
            list_ft.remove(feat_n)
            ttk_aw = ttk_n
        return QgsGeometry.fromPolyline(list_hsl)
