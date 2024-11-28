// ==UserScript==
// @name         PubMed 页面显示期刊影响因子
// @namespace    http://tampermonkey.net/
// @version      0.4
// @description  在 PubMed 和 PMC 页面显示期刊影响因子
// @author       Your name
// @match        https://pubmed.ncbi.nlm.nih.gov/*
// @match        https://pmc.ncbi.nlm.nih.gov/articles/*
// @updateURL    https://raw.githubusercontent.com/ydx1013/DATA/refs/heads/master/pubmed.js
// @downloadURL  https://raw.githubusercontent.com/ydx1013/DATA/refs/heads/master/pubmed.js
// @grant        none
// ==/UserScript==

(function() {
    'use strict';

    // 预定义颜色池
    const colorSchemes = [
        { bg: '#fce4ec', border: '#fce4ec', text: '#c2185b' },  // 粉红
        { bg: '#e3f2fd', border: '#e3f2fd', text: '#1565C0' },  // 蓝色
        { bg: '#e0f7fa', border: '#e0f7fa', text: '#00838F' },  // 青色
        { bg: '#e8f5e9', border: '#e8f5e9', text: '#2E7D32' },  // 绿色
        { bg: '#f3e5f5', border: '#f3e5f5', text: '#6A1B9A' },  // 紫色
        { bg: '#fff3e0', border: '#fff3e0', text: '#E65100' },  // 橙色
        { bg: '#e8eaf6', border: '#e8eaf6', text: '#283593' },  // 靛蓝
        { bg: '#f1f8e9', border: '#f1f8e9', text: '#558b2f' }   // 浅绿
    ];

    // 定义需要显示的指标
    const desiredMetrics = {
        '中科院预警': 'sciwarn',
        '中科院Top分区': 'sciUpTop',
        '中科院基础版': 'sciBase',
        '中科院升级版': 'sciUp',
        '影响因子': 'sciif',
        'JCR分区': 'sci',
        '5年影响因子': 'sciif5',
        'JCI指数': 'jci',
        'ESI学科分类': 'esi'
    };

    // 移除中文字符的函数
    function removeChineseCharacters(str) {
        return str.replace(/[\u4e00-\u9fa5]+/g, '').trim();
    }

    // 获取指标的className
    function getMetricClassName(metric) {
        return metric
            .toLowerCase()
            .replace(/[0-9]/g, '')
            .replace(/最新/g, '')
            .replace(/[：:]/g, '')
            .replace(/\s+/g, '-')
            .replace(/[年]/g, '')
            .trim();
    }

    // 添加缓存系统
    const journalCache = new Map();
    const requestQueue = [];
    let isProcessing = false;
    const RATE_LIMIT = 100; // 每次请求间隔500毫秒


    // 处理请求队列
    async function processQueue() {
        if (isProcessing || requestQueue.length === 0) return;

        isProcessing = true;
        while (requestQueue.length > 0) {
            const { journalName, resolve, reject } = requestQueue.shift();
            try {
                const result = await executeRequest(journalName);
                resolve(result);
            } catch (error) {
                reject(error);
            }
            await new Promise(resolve => setTimeout(resolve, RATE_LIMIT));
        }
        isProcessing = false;
    }

    // 执行实际的API请求
async function executeRequest(journalName) {
    try {
        // 移除URL末尾的斜杠
        const url = 'https://journal-api.yueyang.eu.org';  // 修改这里
        const response = await fetch(url, {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({
                journal: journalName
            })
        });

        const data = await response.json();

        if (data.code === 200 && data.data && data.data.officialRank && data.data.officialRank.all) {
            return processJournalData(data.data.officialRank.all);
        }
        return null;
    } catch (error) {
        console.error('获取期刊数据失败:', error);
        return null;
    }
}

    // 获取期刊数据
    async function getJournalData(journalName) {
        // 检查缓存
        if (journalCache.has(journalName)) {
            return journalCache.get(journalName);
        }

        return new Promise((resolve, reject) => {
            requestQueue.push({
                journalName,
                resolve: (data) => {
                    if (data) {
                        journalCache.set(journalName, data);
                    }
                    resolve(data);
                },
                reject
            });
            processQueue();
        });
    }

    // 处理期刊数据
    function processJournalData(metrics) {
        const result = {};
        for (const [displayName, key] of Object.entries(desiredMetrics)) {
            if (metrics[key] && metrics[key].toString().trim() !== '') {
                result[displayName] = metrics[key];
            }
        }
        return Object.keys(result).length > 0 ? result : null;
    }

    // 更新搜索页面的指标显示
    function updateSearchMetricsDisplay(container, journalData) {
        // 添加搜索页面的样式
        const style = document.createElement('style');
        style.textContent = `
            .journal-metrics {
                display: flex;
                flex-wrap: wrap;
                gap: 4px;
                margin-top: 4px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            }
            .metric-item {
                padding: 2px 6px;
                border-radius: 12px;
                font-size: 12px;
                font-weight: 600;  // 添加这行来设置字体加粗

                width: fit-content;
                display: flex;
                align-items: center;
            }
        `;

        // 为每个指标添加样式
        Object.keys(desiredMetrics).forEach((metric, index) => {
            const className = getMetricClassName(metric);
            const colorScheme = colorSchemes[index % colorSchemes.length];
            style.textContent += `
                .metric-item.${className} {
                    background-color: ${colorScheme.bg};
                    border: none;
                    color: ${colorScheme.text};
                }
            `;
        });

        document.head.appendChild(style);

        const metricsContainer = document.createElement('div');
        metricsContainer.className = 'journal-metrics';

        Object.entries(journalData).forEach(([metric, value]) => {
            const metricElement = document.createElement('div');
            const className = getMetricClassName(metric);
            metricElement.className = `metric-item ${className}`;
            metricElement.textContent = `${metric}: ${value}`;
            metricsContainer.appendChild(metricElement);
        });

        container.appendChild(metricsContainer);
    }

    // 更新详情页面的指标显示
    function updateDetailMetricsDisplay(container, journalData) {
        // 添加详情页面的样式
        const style = document.createElement('style');
        style.textContent = `
            .journal-metrics {
                display: flex;
                flex-wrap: wrap;
                gap: 8px;
                margin-top: 8px;
                font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
            }
            .metric-item {
                padding: 3px 6px;
                border-radius: 8px;
                font-size: 0.5em;
                width: fit-content;
                display: flex;
                align-items: center;
            }
        `;

        // 为每个指标添加样式
        Object.keys(desiredMetrics).forEach((metric, index) => {
            const className = getMetricClassName(metric);
            const colorScheme = colorSchemes[index % colorSchemes.length];
            style.textContent += `
                .metric-item.${className} {
                    background-color: ${colorScheme.bg};
                    border: none;
                    color: ${colorScheme.text};
                }
            `;
        });

        document.head.appendChild(style);

        const metricsContainer = document.createElement('div');
        metricsContainer.className = 'journal-metrics';

        Object.entries(journalData).forEach(([metric, value]) => {
            const metricElement = document.createElement('div');
            const className = getMetricClassName(metric);
            metricElement.className = `metric-item ${className}`;
            metricElement.textContent = `${metric}: ${value}`;
            metricsContainer.appendChild(metricElement);
        });

        container.appendChild(metricsContainer);
    }

    // 修改搜索页面的显示逻辑
    async function addImpactFactorsToSearch() {
        const articles = document.querySelectorAll('.docsum-content');
        const promises = [];

        for (const article of articles) {
            const citationSpan = article.querySelector('.docsum-journal-citation');
            if (!citationSpan ||
                citationSpan.querySelector('.journal-metrics') ||
                citationSpan.dataset.processing === 'true') continue;

            const citationText = citationSpan.textContent;
            const journalRegex = /^([^.]+)\./;
            const match = citationText.match(journalRegex);

            if (match && match[1]) {
                const journalName = removeChineseCharacters(match[1]);
                citationSpan.dataset.processing = 'true';

                promises.push(
                    getJournalData(journalName).then(journalData => {
                        if (journalData) {
                            updateSearchMetricsDisplay(citationSpan, journalData);
                        }
                    }).catch(error => {
                        console.error('处理期刊数据失败:', error);
                    }).finally(() => {
                        delete citationSpan.dataset.processing;
                    })
                );
            }
        }

        await Promise.all(promises);
    }

    // 修改详情页面的显示逻辑
    async function addImpactFactorsToDetail() {
        const journalElement = document.getElementById('full-view-journal-trigger');
        const existingMetrics = document.querySelector('.journal-metrics');

        if (!journalElement ||
            existingMetrics ||
            journalElement.dataset.processing === 'true') return;

        const journalName = removeChineseCharacters(journalElement.textContent);
        journalElement.dataset.processing = 'true';

        try {
            const journalData = await getJournalData(journalName);
            if (journalData) {
                const titleElement = document.querySelector('.heading-title');
                if (titleElement) {
                    updateDetailMetricsDisplay(titleElement, journalData);
                }
            }
        } catch (error) {
            console.error('处理期刊数据失败:', error);
        } finally {
            delete journalElement.dataset.processing;
        }
    }

    // 添加PMC页面的期刊信息处理
    async function addImpactFactorsToPMC() {
        const journalElement = document.querySelector('main article section:first-child div button');
        const titleContainer = document.querySelector('main article section:nth-child(2) div hgroup');

        if (!journalElement || !titleContainer || titleContainer.querySelector('.journal-metrics')) {
            return;
        }


        const journalName = removeChineseCharacters(journalElement.textContent);

        try {
            const journalData = await getJournalData(journalName);
            if (journalData) {
                // 添加PMC页面的样式
                const style = document.createElement('style');
                style.textContent = `
                    .journal-metrics {
                        display: flex;
                        flex-wrap: wrap;
                        gap: 8px;
                        margin-top: 8px;
                        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Arial, sans-serif;
                    }
                    .metric-item {
                        padding: 3px 8px;
                        border-radius: 8px;
                        font-size: 0.9em;
                        width: fit-content;
                        display: flex;
                        align-items: center;
                    }
                `;

                // 为每个指标添加样式
                Object.keys(desiredMetrics).forEach((metric, index) => {
                    const className = getMetricClassName(metric);
                    const colorScheme = colorSchemes[index % colorSchemes.length];
                    style.textContent += `
                        .metric-item.${className} {
                            background-color: ${colorScheme.bg};
                            border: none;
                            color: ${colorScheme.text};
                        }
                    `;
                });

                document.head.appendChild(style);

                const metricsContainer = document.createElement('div');
                metricsContainer.className = 'journal-metrics';

                Object.entries(journalData).forEach(([metric, value]) => {
                    const metricElement = document.createElement('div');
                    const className = getMetricClassName(metric);
                    metricElement.className = `metric-item ${className}`;
                    metricElement.textContent = `${metric}: ${value}`;
                    metricsContainer.appendChild(metricElement);
                });

                titleContainer.appendChild(metricsContainer);
            }
        } catch (error) {
            console.error('处理PMC期刊数据失败:', error);
        }
    }

    // 主函数
    async function init() {
        const isPMC = window.location.hostname.includes('pmc.ncbi.nlm.nih.gov');
        const isPubMedDetail = window.location.pathname.match(/\/\d+\//);

        if (isPMC) {
            await addImpactFactorsToPMC();
        } else if (isPubMedDetail) {
            await addImpactFactorsToDetail();
        } else {
            await addImpactFactorsToSearch();
        }
    }

    // 在页面加载完成后执行
    document.addEventListener('DOMContentLoaded', init);

    // 创建观察器处理动态加载的内容
    const observer = new MutationObserver((mutations) => {
        const isPMC = window.location.hostname.includes('pmc.ncbi.nlm.nih.gov');
        if (isPMC) {
            const titleContainer = document.querySelector('main article section:nth-child(2) div hgroup');
            if (titleContainer && !titleContainer.querySelector('.journal-metrics')) {
                setTimeout(init, 100);
            }
        } else {
            setTimeout(init, 100);
        }
    });

    // 观察页面变化
    const target = document.querySelector('.search-results-chunk') ||
                  document.querySelector('main') ||
                  document.body;

    if (target) {
        observer.observe(target, {
            childList: true,
            subtree: true
        });
    }
})();
